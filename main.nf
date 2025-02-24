nextflow.enable.dsl=2

include { printHeader; helpMessage } from './help' params ( params )
include { RNASEQ_WF } from './workflow/rnaseq.nf' params ( params )


if ( params.help ) {
    helpMessage()
    exit 0
}



params.salmon_index               = WorkflowscRNASeq.getGenomeAttribute( params, 'salmon_index' )
params.tx2gene                    = WorkflowscRNASeq.getGenomeAttribute( params, 'tx2gene' )
params.star_index                 = WorkflowscRNASeq.getGenomeAttribute( params, 'star_index' )
params.gtf_file                   = WorkflowscRNASeq.getGenomeAttribute( params, 'gtf_file' )
params.genebody_ref               = params.genomes [ params.genome ] [ 'genebody_ref']

workflow {
    
    printHeader()
    if ( params.reads ) {
                
        ch_reads = Channel.fromFilePairs( params.reads , size: -1 , checkExists: true )
        ch_reads.ifEmpty{ exit 1, "ERROR: cannot find any fastq files matching the pattern: ${params.reads}\nMake sure that the input file exists!" }

    } else if ( params.input_csv  ) {
                
        ch_reads = Channel.fromPath( params.input_csv  ).splitCsv( header:true )
                        .map { row -> [ row.biosampleName, [ row.read1, row.read2 ] ] }
        ch_reads.ifEmpty{ exit 1, "ERROR: Input csv file is empty." }
    }
                        
     ch_multiqc_config = Channel.fromPath(params.multiqc_config, checkIfExists: true)
     ch_reference_celltype = Channel.fromPath(params.celltype_ref)    

    ch_input_csv = file(params.input_csv, checkIfExists: true)
       
    RNASEQ_WF( 
                params.publish_dir,
                params.enable_publish,
                params.disable_publish,
                ch_reads, 
                ch_input_csv,
                params.adapter_sequence,
                params.adapter_sequence_r2,
                params.salmon_index,
                params.star_index, 
                params.gtf_file, 
                params.tx2gene,
                ch_reference_celltype,
                ch_multiqc_config,
                params.project,
                params.genebody_ref,
                params.min_reads

             )
}

// OnComplete
workflow.onComplete{
    println( "\nPipeline completed successfully.\n\n" )
}

// OnError
workflow.onError{
    println( "\nPipeline failed.\n\n" )
}