params.parameters = ''
params.njobs = 10
params.outdir = 'output'

// Process parameters
parameters = file(params.parameters)

// Read parameters
// PARAMS = Channel.fromPath(paramters).
//     splitCsv(header:false).
//     map{row -> tuple(row[0], row[1], row[2])}
PARAMS = parameters.readLines()

process ED{
    publishDir params.output
    maxForks params.njobs

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
