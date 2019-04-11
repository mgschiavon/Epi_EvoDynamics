// nextflow pipeline for ED_run.out
// </PATH/TO/>nextflow run </PATH/TO/>ED_run.nf --parameters <myparams_file>
//                                      --njonbs <#jobs> --outdir <output_dir>
//                                      -with-trace -with-timeline -with-report &> log &

// Command line options
params.parameters = ''
params.njobs = 10
params.outdir = 'output'

// Process parameters
parameters = file(params.parameters)

// Read parameters
PARAMS = parameters.readLines()

process ED{
    publishDir params.outdir
    maxForks params.njobs
    errorStrategy 'ignore'

    input:
    val params from PARAMS

    output:
    file "*.dat"

    """
    ${workflow.projectDir}/ED_run.out $params
    """
}

// nextflow.config template
/*
process{
    executor = 'local'
    // executor = 'slurm'
}
*/
