#!/usr/bin/env nextflow
import groovy.json.JsonOutput
nextflow.enable.dsl = 2

// Every dependency below is a staged path: publication is never used as a channel.
def stageCommand(stage, key, output, inputs, dependencies, scripts, helpers, fixture, cpus, memoryMb, manifest = '', sourceHash = '') {
    def paths = { value -> (value instanceof Collection ? value : [value]).collect { it.toString() } }
    def configuration = [stage:stage, key:key, output:output, inputs:paths(inputs),
        dependencies:paths(dependencies), scripts:scripts.toString(), fixture:fixture.toString(),
        manifest:manifest.toString(), cpus:cpus, memory_mb:memoryMb,
        engine:params.engine, container_image:params.container_image,
        container_platform:params.container_platform, asmcode:params.asmcode,
        threads:params.threads, memory:params.memory, correction_de:params.correction_de,
        pvalue:params.pvalue, operons_detection:(stage in ['rna_analyze', 'integrate'] ? params.operons_detection : false),
        operons_threshold:params.operons_threshold, amino_diff:params.amino_diff,
        similarity:params.similarity, multitaxa:params.multitaxa,
        max_gse_parallel:params.max_gse_parallel, export_tsv:params.export_tsv,
        mock_delay:params.mock_delay,
        mock_fail_stage:(stage == params.mock_fail_stage && (!params.mock_fail_gse || key == params.mock_fail_gse) ? stage : ''),
        mock_fail_gse:(stage == params.mock_fail_stage && key == params.mock_fail_gse ? key : '')]
    if (sourceHash) configuration.source_sha256 = sourceHash
    """
    cat > task.json <<'PROMOTHEUS_TASK_CONFIG'
    ${JsonOutput.toJson(configuration)}
    PROMOTHEUS_TASK_CONFIG
    python3 ${helpers}/run_modular_stage.py --config task.json
    """
}

process VALIDATE_INPUTS {
    tag 'validate and split GSE inputs'
    publishDir "${params.output_dir}/.promotheus-nextflow", mode:'copy', pattern:'input_*', overwrite:true
    input:
    path input_files
    path helpers
    val source_hash
    output:
    tuple path('input_manifest.tsv'), path('input_validation.json'), emit: validated
    path 'shards', emit: shards
    script:
    """
    # Source content SHA-256: ${source_hash}
    mkdir source_input
    cp *.xlsx source_input/
    python3 ${helpers}/validate_promotheus_inputs.py --input-dir source_input \
      --asmcode '${params.asmcode}' --manifest input_manifest.tsv --report input_validation.json
    python3 ${helpers}/prepare_gse_inputs.py --input-dir source_input --output-dir shards
    """
}

// Conservative dependency coverage: include every adapter/scientific resource,
// including nested imports and the Rockhopper JAR; ignore generated cache files.
def sourceCodeHash(root, folders) {
    def rootPath = new File(root.toString()).toPath()
    def sources = []
    folders.each { folder ->
        new File(root.toString(), folder).eachFileRecurse { file ->
            def relative = rootPath.relativize(file.toPath()).toString()
            if (file.isFile() && !relative.tokenize('/').any { it.startsWith('.') || it == '__pycache__' } && !relative.endsWith('.pyc'))
                sources.add(relative)
        }
    }
    def hash = java.security.MessageDigest.getInstance('SHA-256')
    sources.sort().each { name ->
        hash.update(name.getBytes('UTF-8'))
        hash.update((byte)0)
        hash.update(java.security.MessageDigest.getInstance('SHA-256').digest(new File(root.toString(), name).bytes))
    }
    hash.digest().encodeHex().toString()
}

def orthologyForkLimit() {
    def requested = params.max_orthology_parallel as int
    if (!params.reserve_rnaseq_resources) return requested
    if ((params.max_gse_parallel as int) < 1 || (params.orthology_cpus as int) < 1 || (params.orthology_memory as int) < 512)
        error 'Balanced scheduling requires positive concurrency/CPU requests and at least 512 MB per orthology pair'
    def taskCpus = params.task_cpus ?: Math.max(1, (int)((params.threads as int) / (params.max_gse_parallel as int)))
    def taskMemory = params.task_memory ?: Math.max(512, (int)((params.memory as int) / (params.max_gse_parallel as int)))
    def cpuSlots = (int)Math.floor(((params.threads as int) - (taskCpus as int)) / (double)(params.orthology_cpus as int))
    def memorySlots = (int)Math.floor(((params.memory as int) - (taskMemory as int)) / (double)(params.orthology_memory as int))
    def effective = Math.min(requested, Math.min(cpuSlots, memorySlots))
    if (effective < 1)
        error 'Cannot reserve one RNA-Seq task alongside one orthology pair: increase --threads/--memory or reduce per-task requests; disable --reserve_rnaseq_resources for orthology-only runs'
    return effective
}

def validateOrthologyResources() {
    if ((params.max_orthology_parallel as int) < 1) error '--max_orthology_parallel must be positive'
    if ((params.orthology_cpus as int) < 1 || (params.orthology_cpus as int) > (params.threads as int))
        error '--orthology_cpus must be positive and fit within --threads'
    if ((params.orthology_memory as int) < 512 || (params.orthology_memory as int) > (params.memory as int))
        error '--orthology_memory must be at least 512 MB and fit within --memory'
    def forks = orthologyForkLimit()
    if (params.reserve_rnaseq_resources)
        log.info "RNA-Seq resource reserve enabled: orthology concurrency capped at ${forks} (requested ${params.max_orthology_parallel})"
}

def configuredOrthologyForks = orthologyForkLimit()

process ORTHOLOGY_PREPARE {
    label 'scientific'
    tag 'download and prepare reference proteomes'
    input:
    path input_files
    tuple path(manifest), path(validation)
    path genome_cache
    path scripts
    path helpers
    path fixture
    val source_hash
    output:
    path 'reference'
    script:
    stageCommand('orthology_prepare', 'all', 'reference', input_files, genome_cache, scripts, helpers, fixture, task.cpus, task.memory.toMega(), '', source_hash)
}

process ORTHOLOGY_PAIR {
    label 'orthology_pair'
    tag "${key}"
    maxForks configuredOrthologyForks
    input:
    tuple val(key), path(pair)
    path shared
    path scripts
    path helpers
    path fixture
    val source_hash
    output:
    path "pair_${key}"
    script:
    stageCommand('orthology_pair', key, "pair_${key}", pair, shared, scripts, helpers, fixture, task.cpus, task.memory.toMega(), '', source_hash)
}

process ORTHOLOGY_MERGE {
    label 'scientific'
    tag 'merge orthology pairs in reference order'
    publishDir "${params.output_dir}/orthology", mode:'copy', pattern:'orthologous_genes.txt'
    input:
    path reference
    path mappings
    path scripts
    path helpers
    path fixture
    val source_hash
    output:
    path 'orthology', emit: bundle
    path 'orthologous_genes.txt', emit: table
    script:
    stageCommand('orthology_merge', 'all', 'orthology', [], [reference] + mappings, scripts, helpers, fixture, task.cpus, task.memory.toMega(), '', source_hash) + """
    cp orthology/orthologous_genes.txt orthologous_genes.txt
    """
}

workflow ORTHOLOGY {
    take:
    inputs
    validated
    genome_cache
    scripts
    helpers
    fixture
    source_hash
    main:
    ORTHOLOGY_PREPARE(inputs, validated, genome_cache, scripts, helpers, fixture, source_hash)
    pairs = ORTHOLOGY_PREPARE.out.flatMap { root ->
        root.resolve('pairs.tsv').readLines().drop(1).collect { line ->
            def fields = line.split('\t')
            tuple(fields[1], root.resolve("pairs/${fields[1]}"))
        }
    }
    shared = ORTHOLOGY_PREPARE.out.map { it.resolve('shared') }
    ORTHOLOGY_PAIR(pairs, shared, scripts, helpers, fixture, source_hash)
    ORTHOLOGY_MERGE(ORTHOLOGY_PREPARE.out, ORTHOLOGY_PAIR.out.toList(), scripts, helpers, fixture, source_hash)
    emit:
    result = ORTHOLOGY_MERGE.out.bundle
}

workflow ORTHOLOGY_ONLY {
    validateOrthologyResources()
    inputs = Channel.fromPath("${params.input_dir}/*.xlsx", checkIfExists:true).collect()
    helpers = file("${projectDir}/bin", checkIfExists:true)
    scripts = file("${projectDir}/Scripts", checkIfExists:true)
    fixture = file("${projectDir}/tests/fixtures/toy_expected/final_table.xlsx", checkIfExists:true)
    genome_cache = params.orthology_genomes ? file(params.orthology_genomes, checkIfExists:true) : []
    helper_hash = sourceCodeHash(projectDir, ['bin'])
    scientific_hash = sourceCodeHash(projectDir, ['bin', 'Scripts'])
    VALIDATE_INPUTS(inputs, helpers, helper_hash)
    ORTHOLOGY(inputs, VALIDATE_INPUTS.out.validated, genome_cache, scripts, helpers, fixture, scientific_hash)
}

process MICROARRAY_GSE {
    label 'scientific'
    tag "${accession}"
    maxForks params.max_gse_parallel
    input:
    tuple val(accession), val(group), val(ordinal), path(shard)
    path orthology
    path scripts
    path helpers
    path fixture
    val source_hash
    output:
    path "microarray_${accession}"
    script:
    stageCommand('microarray', accession, "microarray_${accession}", shard, orthology, scripts, helpers, fixture, task.cpus, task.memory.toMega(), '', source_hash)
}

process RNA_DOWNLOAD_GSE {
    label 'scientific'
    tag "${accession}"
    maxForks params.max_gse_parallel
    input:
    tuple val(accession), val(group), val(ordinal), path(shard)
    path scripts
    path helpers
    path fixture
    val source_hash
    output:
    tuple val(accession), val(group), val(ordinal), path("download_${accession}")
    script:
    stageCommand('rna_download', accession, "download_${accession}", shard, [], scripts, helpers, fixture, task.cpus, task.memory.toMega(), '', source_hash)
}

process RNA_ALIGN_GSE {
    label 'scientific'
    tag "${accession}"
    maxForks params.max_gse_parallel
    input:
    tuple val(accession), val(group), val(ordinal), path(download)
    path scripts
    path helpers
    path fixture
    val source_hash
    output:
    tuple val(accession), val(group), val(ordinal), path("aligned_${accession}")
    script:
    stageCommand('rna_align', accession, "aligned_${accession}", [], download, scripts, helpers, fixture, task.cpus, task.memory.toMega(), '', source_hash)
}

process RNA_ANALYZE {
    label 'scientific'
    tag "${key}"
    maxForks params.max_gse_parallel
    input:
    tuple val(key), path(aligned)
    path scripts
    path helpers
    path fixture
    val source_hash
    output:
    path "rna_${key}"
    script:
    stageCommand('rna_analyze', key, "rna_${key}", [], aligned, scripts, helpers, fixture, task.cpus, task.memory.toMega(), '', source_hash)
}

process INTEGRATE_RESULTS {
    label 'scientific'
    tag 'integrate all GSE results'
    publishDir "${params.output_dir}", mode:'copy', pattern:'{results,RNA_Seq,orthologous_genes.txt}', overwrite:true
    publishDir "${params.output_dir}/.promotheus-nextflow", mode:'copy', pattern:'{run_record.tsv,input_manifest.tsv,input_validation.json,checkpoint_manifest.tsv}', overwrite:true
    input:
    path input_files
    tuple path(manifest), path(validation)
    path orthology
    path microarray
    path rnaseq
    path scripts
    path helpers
    path fixture
    val source_hash
    output:
    tuple path('results'), path('RNA_Seq'), path('orthologous_genes.txt'), path('run_record.tsv'), path('input_manifest.tsv'), path('input_validation.json'), path('checkpoint_manifest.tsv')
    script:
    stageCommand('integrate', 'all', 'integrated', input_files, [orthology] + microarray + rnaseq, scripts, helpers, fixture, task.cpus, task.memory.toMega(), manifest, source_hash) + """
    mv integrated/results integrated/RNA_Seq integrated/orthologous_genes.txt .
    mv integrated/.promotheus-nextflow/run_record.tsv integrated/.promotheus-nextflow/checkpoint_manifest.tsv .
    cp -L ${manifest} manifest_copy.tsv
    cp -L ${validation} validation_copy.json
    mv -f manifest_copy.tsv input_manifest.tsv
    mv -f validation_copy.json input_validation.json
    """
}

process VERIFY_OUTPUTS {
    tag 'verify integrated results'
    publishDir "${params.output_dir}/.promotheus-nextflow", mode:'copy', overwrite:true
    input:
    tuple path(results), path(rnaseq), path(orthologues), path(record), path(manifest), path(validation), path(checkpoints)
    path helpers
    val source_hash
    output:
    tuple path('output_manifest.tsv'), path('output_validation.json')
    script:
    """
    # Source content SHA-256: ${source_hash}
    python3 ${helpers}/verify_promotheus_outputs.py --output-dir . \
      --run-record ${record} \
      --input-manifest ${manifest} \
      --checkpoint-manifest ${checkpoints} \
      --manifest output_manifest.tsv --report output_validation.json
    """
}

workflow {
    validateOrthologyResources()
    if (!params.input_dir || !params.asmcode) error 'Required parameters: --input_dir and --asmcode'
    if (!(params.asmcode ==~ /ASM\d+v\d+/)) error '--asmcode must match ASM<digits>v<digits>'
    if (!(params.engine in ['mock', 'docker'])) error '--engine must be docker or mock'
    if ((params.threads as int) < 1 || (params.max_gse_parallel as int) < 1) error 'threads and max_gse_parallel must be positive'
    if ((params.memory as int) < 512) error '--memory must be at least 512 MB'
    if (params.task_cpus && (params.task_cpus as int) > (params.threads as int)) error '--task_cpus exceeds --threads'
    if (params.task_memory && (params.task_memory as int) > (params.memory as int)) error '--task_memory exceeds --memory'
    if (params.task_cpus != null && (params.task_cpus as int) < 1) error '--task_cpus must be positive'
    if (params.task_memory != null && (params.task_memory as int) < 512) error '--task_memory must be at least 512 MB'
    if (params.cache_dir) error 'The modular workflow uses Nextflow work directories and -resume; --cache_dir is supported by main_legacy.nf only'
    inputs = Channel.fromPath("${params.input_dir}/*.xlsx", checkIfExists:true).collect()
    helpers = file("${projectDir}/bin", checkIfExists:true)
    scripts = file("${projectDir}/Scripts", checkIfExists:true)
    fixture = file("${projectDir}/tests/fixtures/toy_expected/final_table.xlsx", checkIfExists:true)
    helper_hash = sourceCodeHash(projectDir, ['bin'])
    scientific_hash = sourceCodeHash(projectDir, ['bin', 'Scripts'])
    VALIDATE_INPUTS(inputs, helpers, helper_hash)
    datasets = VALIDATE_INPUTS.out.shards.flatMap { root ->
        root.resolve('datasets.tsv').readLines().drop(1).collect { line ->
            def fields = line.split('\t')
            tuple(fields[0], fields[1], fields[2], fields[3] as int, root.resolve(fields[1]))
        }
    }.branch {
        microarray: it[0] == 'microarray'
        rnaseq: it[0] == 'rnaseq'
    }
    ma = datasets.microarray.map { kind, accession, group, ordinal, shard -> tuple(accession, group, ordinal, shard) }
    rna = datasets.rnaseq.map { kind, accession, group, ordinal, shard -> tuple(accession, group, ordinal, shard) }
    genome_cache = params.orthology_genomes ? file(params.orthology_genomes, checkIfExists:true) : []
    ORTHOLOGY(inputs, VALIDATE_INPUTS.out.validated, genome_cache, scripts, helpers, fixture, scientific_hash)
    MICROARRAY_GSE(ma, ORTHOLOGY.out, scripts, helpers, fixture, scientific_hash)
    RNA_DOWNLOAD_GSE(rna, scripts, helpers, fixture, scientific_hash)
    RNA_ALIGN_GSE(RNA_DOWNLOAD_GSE.out, scripts, helpers, fixture, scientific_hash)
    if (params.operons_detection) {
        analysis = RNA_ALIGN_GSE.out.map { accession, group, ordinal, aligned -> tuple(group, aligned) }.groupTuple()
    } else {
        analysis = RNA_ALIGN_GSE.out.map { accession, group, ordinal, aligned -> tuple(accession, [aligned]) }
    }
    RNA_ANALYZE(analysis, scripts, helpers, fixture, scientific_hash)
    INTEGRATE_RESULTS(inputs, VALIDATE_INPUTS.out.validated, ORTHOLOGY.out,
                      MICROARRAY_GSE.out.collect(), RNA_ANALYZE.out.collect(), scripts, helpers, fixture, scientific_hash)
    VERIFY_OUTPUTS(INTEGRATE_RESULTS.out, helpers, helper_hash)
}
