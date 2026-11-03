#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input_dir          = params.input_dir ?: null
params.output_dir         = params.output_dir ?: 'results'
params.container_image    = params.container_image ?: 'ste40/promotheus@sha256:a7c5196b4ebf6f93ef8a4ff2698481d24e3f63306c9eaccd413bce3c48eebec7'
params.container_platform = params.container_platform ?: 'linux/amd64'
params.engine             = params.engine ?: 'docker'
params.threads            = params.threads ?: 8
params.memory             = params.memory ?: 4000
params.correction_de      = params.correction_de ?: 'none'
params.pvalue             = params.pvalue ?: 0.05
params.operons_detection  = params.operons_detection == null ? true : params.operons_detection
params.amino_diff         = params.amino_diff ?: 3
params.similarity         = params.similarity ?: 30
params.asmcode            = params.asmcode ?: null
params.operons_threshold  = params.operons_threshold ?: 3
params.multitaxa          = params.multitaxa == null ? false : params.multitaxa
params.max_gse_parallel   = params.max_gse_parallel ?: 2
params.cache_dir          = params.cache_dir ?: null
params.keep_cache         = params.keep_cache == null ? true : params.keep_cache
params.export_tsv         = params.export_tsv == null ? true : params.export_tsv

process VALIDATE_INPUTS {
    tag 'validate Promotheus inputs'

    publishDir "${params.output_dir}/.promotheus-nextflow", mode: 'copy', overwrite: true

    input:
    path input_files
    path validator

    output:
    tuple path('input_manifest.tsv'), path('input_validation.json'), emit: validated

    script:
    def asmcode = params.asmcode ?: ''
    """
    python3 ${validator} \
      --input-dir . \
      --asmcode '${asmcode}' \
      --manifest input_manifest.tsv \
      --report input_validation.json
    """
}

process RUN_PROMOTHEUS {
    tag "Promotheus (${params.engine})"

    cpus { params.threads as int }
    memory { "${params.memory} MB" }
    publishDir "${params.output_dir}/.promotheus-nextflow", mode: 'copy', overwrite: true

    input:
    tuple path(input_manifest), path(input_validation)
    path runner
    path toy_final_table
    path orthology_script
    path refseq_download_script
    path microarray_metadata_script
    path microarray_intra_script
    path rna_data_retrieval_script
    path rna_bam_parsing_script
    path rna_nde_detection_script
    path rna_pipeline_script
    path rna_operons_script
    path interexp_analysis_script
    path xlsx_repairer
    path xlsx_integrity_library
    path tsv_exporter
    path checkpoint_builder

    output:
    tuple path('run_record.tsv'), path('input_manifest.tsv'), path('checkpoint_manifest.tsv'), emit: record

    script:
    def inputDir = file(params.input_dir).toAbsolutePath()
    def outputDir = file(params.output_dir).toAbsolutePath()
    def cacheDir = params.cache_dir \
        ? file(params.cache_dir).toAbsolutePath() \
        : file("${params.output_dir}/.promotheus-cache").toAbsolutePath()
    def fixtureArgument = params.engine == 'mock' ? "--fixture-output '${toy_final_table}'" : ''
    def refseqPatchArguments = params.engine == 'docker' ? """\
      --orthology-script '${orthology_script}' \\
      --refseq-download-script '${refseq_download_script}' \\
      --microarray-metadata-script '${microarray_metadata_script}' \\
      --microarray-intra-script '${microarray_intra_script}' \\
      --rna-data-retrieval-script '${rna_data_retrieval_script}' \\
      --rna-bam-parsing-script '${rna_bam_parsing_script}' \\
      --rna-nde-detection-script '${rna_nde_detection_script}' \\
      --rna-pipeline-script '${rna_pipeline_script}' \\
      --rna-operons-script '${rna_operons_script}' \\
      --interexp-analysis-script '${interexp_analysis_script}'""" : ''
    """
    bash ${runner} \
      --input-dir '${inputDir}' \
      --output-dir '${outputDir}' \
      --input-manifest '${input_manifest}' \
      --record run_record.tsv \
      --image '${params.container_image}' \
      --platform '${params.container_platform}' \
      --threads '${params.threads}' \
      --memory '${params.memory}' \
      --correction-de '${params.correction_de}' \
      --pvalue '${params.pvalue}' \
      --operons-detection '${params.operons_detection}' \
      --amino-diff '${params.amino_diff}' \
      --similarity '${params.similarity}' \
      --asmcode '${params.asmcode}' \
      --operons-threshold '${params.operons_threshold}' \
      --multitaxa '${params.multitaxa}' \
      --max-gse-parallel '${params.max_gse_parallel}' \
      --cache-dir '${cacheDir}' \
      --keep-cache '${params.keep_cache}' \
      --export-tsv '${params.export_tsv}' \
      ${refseqPatchArguments} \
      ${fixtureArgument}
    python3 ${xlsx_repairer} \
      --workbook '${outputDir}/results/final_table.xlsx'
    if [[ '${params.export_tsv}' == 'true' ]]; then
      python3 ${tsv_exporter} \
        --workbook '${outputDir}/results/final_table.xlsx' \
        --output '${outputDir}/results/final_table.tsv'
    fi
    python3 ${checkpoint_builder} \
      --output-dir '${outputDir}' \
      --input-manifest '${input_manifest}' \
      --manifest checkpoint_manifest.tsv
    """
}

process VERIFY_OUTPUTS {
    tag 'verify Promotheus outputs'

    publishDir "${params.output_dir}/.promotheus-nextflow", mode: 'copy', overwrite: true

    input:
    tuple path(run_record), path(input_manifest), path(checkpoint_manifest)
    path verifier
    path validator_library
    path xlsx_integrity_library

    output:
    tuple path('output_manifest.tsv'), path('output_validation.json'), emit: verified

    script:
    def outputDir = file(params.output_dir).toAbsolutePath()
    """
    python3 ${verifier} \
      --output-dir '${outputDir}' \
      --run-record '${run_record}' \
      --input-manifest '${input_manifest}' \
      --checkpoint-manifest '${checkpoint_manifest}' \
      --manifest output_manifest.tsv \
      --report output_validation.json
    """
}

workflow {
    if (!params.input_dir) {
        error 'Missing required parameter: --input_dir'
    }
    if (!params.asmcode) {
        error 'Missing required parameter: --asmcode'
    }
    if (!(params.engine in ['docker', 'mock'])) {
        error "--engine must be 'docker' or 'mock'; received: ${params.engine}"
    }
    if ((params.max_gse_parallel as int) < 1) {
        error "--max_gse_parallel must be at least 1; received: ${params.max_gse_parallel}"
    }

    input_dir = file(params.input_dir, checkIfExists: true)
    if (!input_dir.isDirectory()) {
        error "--input_dir is not a directory: ${input_dir}"
    }

    input_files = Channel
        .fromPath("${input_dir}/*.xlsx", checkIfExists: true)
        .collect()

    validator = file("${projectDir}/bin/validate_promotheus_inputs.py", checkIfExists: true)
    verifier = file("${projectDir}/bin/verify_promotheus_outputs.py", checkIfExists: true)
    runner = params.engine == 'mock' \
        ? file("${projectDir}/tests/bin/run_mock_promotheus.sh", checkIfExists: true) \
        : file("${projectDir}/bin/run_promotheus_docker.sh", checkIfExists: true)
    toy_final_table = file("${projectDir}/tests/fixtures/toy_expected/final_table.xlsx", checkIfExists: true)
    orthology_script = file("${projectDir}/Scripts/Microarray/ortologhi_new.R", checkIfExists: true)
    refseq_download_script = file("${projectDir}/Scripts/Microarray/dowload_file_assembly.R", checkIfExists: true)
    microarray_metadata_script = file("${projectDir}/Scripts/Microarray/metadata_retrieval.R", checkIfExists: true)
    microarray_intra_script = file("${projectDir}/Scripts/Microarray/intra_experiment_analysis.R", checkIfExists: true)
    rna_data_retrieval_script = file("${projectDir}/Scripts/RNA-Seq/Data_retrieval.R", checkIfExists: true)
    rna_bam_parsing_script = file("${projectDir}/Scripts/RNA-Seq/BAM_parsing.R", checkIfExists: true)
    rna_nde_detection_script = file("${projectDir}/Scripts/RNA-Seq/NDE_detection.R", checkIfExists: true)
    rna_pipeline_script = file("${projectDir}/Scripts/RNA-Seq/Pipeline_bash_RNA_Seq.R", checkIfExists: true)
    rna_operons_script = file("${projectDir}/Scripts/RNA-Seq/Operons_detection.R", checkIfExists: true)
    interexp_analysis_script = file("${projectDir}/Scripts/Microarray/interexp_analysis.R", checkIfExists: true)
    xlsx_repairer = file("${projectDir}/bin/repair_promotheus_xlsx.py", checkIfExists: true)
    xlsx_integrity_library = file("${projectDir}/bin/xlsx_integrity.py", checkIfExists: true)
    tsv_exporter = file("${projectDir}/bin/export_promotheus_tsv.py", checkIfExists: true)
    checkpoint_builder = file("${projectDir}/bin/build_checkpoint_manifest.py", checkIfExists: true)

    validated = VALIDATE_INPUTS(input_files, validator)
    run_record = RUN_PROMOTHEUS(
      validated,
      runner,
      toy_final_table,
      orthology_script,
      refseq_download_script,
      microarray_metadata_script,
      microarray_intra_script,
      rna_data_retrieval_script,
      rna_bam_parsing_script,
      rna_nde_detection_script,
      rna_pipeline_script,
      rna_operons_script,
      interexp_analysis_script,
      xlsx_repairer,
      xlsx_integrity_library,
      tsv_exporter,
      checkpoint_builder
    )
    VERIFY_OUTPUTS(run_record, verifier, validator, xlsx_integrity_library)
}
