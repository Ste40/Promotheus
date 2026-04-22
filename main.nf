nextflow.enable.dsl=2

params.input_dir          = params.input_dir ?: null
params.output_dir         = params.output_dir ?: 'nf_output'
params.container_image    = params.container_image ?: 'ste40/promotheus:latest'
params.threads            = params.threads ?: 8
params.memory             = params.memory ?: 4000
params.correction_de      = params.correction_de ?: 'none'
params.pvalue             = params.pvalue ?: 0.05
params.operons_detection  = params.operons_detection ?: true
params.amino_diff         = params.amino_diff ?: 3
params.similarity         = params.similarity ?: 30
params.asmcode            = params.asmcode ?: null
params.operons_threshold  = params.operons_threshold ?: 3
params.multitaxa          = params.multitaxa ?: false
params.skip_execution     = params.skip_execution ?: false

process VALIDATE_INPUTS {
    tag 'validate-inputs'

    output:
    path 'input_summary.txt'

    script:
    def inputDir = params.input_dir
    def outputDir = params.output_dir
    def asm = params.asmcode
    """
    set -euo pipefail

    if [[ -z \"${inputDir}\" ]]; then
      echo 'ERROR: --input_dir is required.' >&2
      exit 1
    fi

    if [[ ! -d \"${inputDir}\" ]]; then
      echo "ERROR: input_dir not found: ${inputDir}" >&2
      exit 1
    fi

    if [[ -z \"${asm}\" ]]; then
      echo 'ERROR: --asmcode is required (format ASM<digits>v<digits>).' >&2
      exit 1
    fi

    if [[ ! \"${asm}\" =~ ^ASM[0-9]+v[0-9]+$ ]]; then
      echo "ERROR: invalid asmcode format '${asm}'." >&2
      exit 1
    fi

    micro_count=0
    rna_count=0
    shopt -s nullglob
    for f in \"${inputDir}\"/*; do
      base=$(basename \"$f\")
      [[ \"$base\" == *Microarray* || \"$base\" == *Microarrays* ]] && micro_count=$((micro_count+1))
      [[ \"$base\" == *RNA_Seq* ]] && rna_count=$((rna_count+1))
    done

    mkdir -p \"${outputDir}\"

    {
      echo "input_dir=${inputDir}"
      echo "output_dir=${outputDir}"
      echo "asmcode=${asm}"
      echo "microarray_candidates=${micro_count}"
      echo "rnaseq_candidates=${rna_count}"
    } > input_summary.txt
    """
}

process RUN_PROMOTHEUS_DOCKER {
    tag 'run-promotheus'

    input:
    path summary_file

    output:
    path 'run_metadata.txt'

    when:
    !params.skip_execution

    script:
    def threads = params.threads
    def memory = params.memory
    def correction = params.correction_de
    def pvalue = params.pvalue
    def operons = params.operons_detection
    def amino = params.amino_diff
    def sim = params.similarity
    def asm = params.asmcode
    def opThreshold = params.operons_threshold
    def multitaxa = params.multitaxa
    def inputDir = params.input_dir
    def outputDir = params.output_dir
    def image = params.container_image

    """
    set -euo pipefail

    docker run --rm \
      -e THREADS=${threads} \
      -e MEMORY=${memory} \
      -e CORRECTION_DE=${correction} \
      -e PVALUE=${pvalue} \
      -e OPERONS_DETECTION=${operons} \
      -e AMINO_DIFF=${amino} \
      -e SIMILARITY=${sim} \
      -e ASMcode=${asm} \
      -e OPERONS_THRESHOLD=${opThreshold} \
      -e MULTITAXA=${multitaxa} \
      -v \"${inputDir}:/input\" \
      -v \"${outputDir}:/output\" \
      ${image}

    {
      echo 'Promotheus docker run completed.'
      echo 'Parameters:'
      cat ${summary_file}
      echo "threads=${threads}"
      echo "memory=${memory}"
      echo "correction_de=${correction}"
      echo "pvalue=${pvalue}"
      echo "operons_detection=${operons}"
      echo "amino_diff=${amino}"
      echo "similarity=${sim}"
      echo "operons_threshold=${opThreshold}"
      echo "multitaxa=${multitaxa}"
      echo "container_image=${image}"
    } > run_metadata.txt
    """
}

workflow {
    summary = VALIDATE_INPUTS()
    RUN_PROMOTHEUS_DOCKER(summary)
}
