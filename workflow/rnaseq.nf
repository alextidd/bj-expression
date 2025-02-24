nextflow.enable.dsl=2
import groovy.json.JsonOutput

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/


//
// MODULES
//

include { PUBLISH_INPUT_DATASET_WF } from '../nf-bioskryb-utils/modules/bioskryb/publish_input_dataset/main.nf' addParams(timestamp: params.timestamp)
include { CUSTOM_FASTQ_MERGE_WF } from '../modules/local/custom_fastq_merge/main.nf' addParams(timestamp: params.timestamp)
// include { PETASUITE_DECOMPRESS_WF } from '../nf-bioskryb-utils/modules/petasuite/decompress/main.nf' addParams(timestamp: params.timestamp)
include { SEQTK_WF } from '../nf-bioskryb-utils/modules/seqtk/sample/main.nf' addParams(timestamp: params.timestamp)
include { FastpFull_WF } from '../nf-bioskryb-utils/modules/fastp/main.nf' addParams(timestamp: params.timestamp)
include { SALMONQUANT } from '../nf-bioskryb-utils/modules/salmon/main.nf' addParams(timestamp: params.timestamp)
include { CREATE_TXIMPORT_SALMON_TX_GENE } from '../modules/local/custom_salmonpostprocess/create_tximport_salmon_tx_gene/main.nf' addParams(timestamp: params.timestamp)
include { MERGE_TXIMPORT_SALMON_TX_GENE_WF } from '../modules/local/custom_salmonpostprocess/merge_tximport_salmon_tx_gene/main.nf' addParams(timestamp: params.timestamp)
include { STARALIGN } from '../nf-bioskryb-utils/modules/star/main.nf' addParams(timestamp: params.timestamp)
include { QUALIMAP_BAMRNA } from '../nf-bioskryb-utils/modules/qualimap/rnaqc/main.nf' addParams(timestamp: params.timestamp)
include { STARGETPRIM_WF } from '../modules/local/custom_stargetprim/main.nf' addParams(timestamp: params.timestamp)
include { HTSEQ_COUNTS } from '../nf-bioskryb-utils/modules/htseq/main.nf' addParams(timestamp: params.timestamp)
include { CREATE_HTSEQ_SUMMARY_WF } from '../modules/local/custom_starpostprocess/create_htseq_summary/main.nf' addParams(timestamp: params.timestamp)
include { MERGE_HTSEQ_SUMMARY_WF } from '../modules/local/custom_starpostprocess/merge_htseq_summary/main.nf' addParams(timestamp: params.timestamp)
include { PLOTTER_PCAHEATMAP_HTSEQ_SUMMARY_WF } from '../modules/local/custom_starpostprocess/plotter_pcaheatmapqc_htseq_summary/main.nf' addParams( imestamp: params.timestamp)
include { CELL_TYPING } from '../modules/local/custom_cell_typing/main.nf' addParams(timestamp: params.timestamp)
include { CREATE_QC_REPORT } from '../modules/local/custom_qc_report/main.nf' addParams(timestamp: params.timestamp)
include { MULTIQC_WF } from '../nf-bioskryb-utils/modules/multiqc/main.nf' addParams(timestamp: params.timestamp)
include { REPORT_VERSIONS_WF } from '../nf-bioskryb-utils/modules/bioskryb/report_tool_versions/main.nf' addParams(timestamp: params.timestamp)
// include { CUSTOM_PIPELINE_JSON_OUTPUT } from '../modules/local/custom_pipeline_json_output/main.nf' addParams(timestamp: params.timestamp)
include { CREATE_HTSEQ_MATRIX } from '../modules/local/custom_starpostprocess/create_htseq_matrix/main.nf' addParams(timestamp: params.timestamp)
include { CREATE_MATRIX_SALMON_TX_GENE } from '../modules/local/custom_salmonpostprocess/create_matrix_salmon_tx_gene/main.nf' addParams(timestamp: params.timestamp)
include { GENE_BODY_COVERAGE_RNA } from '../modules/local/gene_body_rna/main.nf' addParams(timestamp: params.timestamp)
include { GENE_BODY_COVERAGE_RNA_PLOT } from '../modules/local/gene_body_rna_plot_all/main.nf' addParams(timestamp: params.timestamp)
include { CREATE_MASTER_STATS } from '../modules/local/create_master_stats/main.nf' addParams(timestamp: params.timestamp)
include { ALEVIN_NOQUANT_DUMPFQ } from  '../nf-bioskryb-utils/modules/alevin/alevin_noquant_dumpfq_10x/main.nf' addParams(timestamp: params.timestamp)
include { CUSTOM_AWK_DEMUX_CBC_FASTQ } from  '../modules/local/custom_10x_demux/awk_demux_cbc_fastq/main.nf' addParams(timestamp: params.timestamp)
include { CALC_DYNAMICRANGE } from '../modules/local/custom_dynamicrange_report/main.nf' addParams(timestamp: params.timestamp)
include { COUNT_READS_FASTQ_WF } from '../nf-bioskryb-utils/modules/bioskryb/custom_read_counts/main.nf' addParams(timestamp: params.timestamp)


//include { STARFUSION_WF } from  '../nf-bioskryb-utils/modules/starFusion/main.nf' addParams(timestamp: params.timestamp)
//include { MERGE_STAR_FUSION } from '../modules/local/custom_starpostprocess/merge_starfusion/main.nf' addParams(timestamp: params.timestamp)




workflow RNASEQ_WF {
    
    take:
        ch_publish_dir
        ch_enable_publish
        ch_disable_publish
        ch_reads
        ch_input_csv 
        ch_adapter_sequence
        ch_adapter_sequence_r2
        ch_salmon_index
        ch_star_index
        ch_gtf
        ch_tx2gene
        ch_reference_celltype
        ch_multiqc_config
        ch_project_name
        ch_genebody_ref
        ch_min_reads

       
    
    main:
    
    if( ! params.skip_10X ){
    
      ALEVIN_NOQUANT_DUMPFQ ( ch_reads, 
                                ch_salmon_index,
                                ch_tx2gene,
                                params.lib_type_10x,
                                params.lib_protocol_10x,
                                ch_publish_dir,
                                ch_disable_publish
       )
       
        CUSTOM_AWK_DEMUX_CBC_FASTQ ( ALEVIN_NOQUANT_DUMPFQ.out.fastq_cbc,
                                ch_publish_dir,
                                ch_disable_publish
        
                              )
                              
        ch_fastqs = CUSTOM_AWK_DEMUX_CBC_FASTQ.out.fastqs
        
        
    }else{
    
        if ( params.input_csv  ) {
            
            PUBLISH_INPUT_DATASET_WF (
                                        ch_input_csv,
                                        ch_publish_dir,
                                        ch_enable_publish
                                    )

        } else {
            ch_input_csv = Channel.of(file("NO_FILE"))
        }
        COUNT_READS_FASTQ_WF (
                                    ch_reads,
                                    params.publish_dir,
                                    params.enable_publish
                                )
        COUNT_READS_FASTQ_WF.out.read_counts
            .map { sample_id, files, read_count -> 
                [sample_id, files, read_count.toInteger()]
            }
            .branch {
                small: it[2] < ch_min_reads
                large: it[2] >= ch_min_reads
            }
            .set { branched_reads }
        
        ch_fastqs = ch_reads
    }
    
    ch_seqtk_version = Channel.empty()
    ch_fastp_report = Channel.empty()
    ch_fastp_version = Channel.empty()
    
    if ( !params.skip_subsampling ) {
        // Run sub sampling followed by fastp only when subsampling is enabled 
        // and when there are reads with min threshold
        if (branched_reads.large.ifEmpty { false }) {
            ch_fastqs_with_nreads = ch_fastqs.map { sample_name, reads ->
                tuple(sample_name, reads, params.n_reads)
            }
            SEQTK_WF (ch_fastqs_with_nreads, false, params.read_length, params.seqtk_sample_seed, ch_publish_dir, ch_disable_publish)
            ch_seqtk_version = SEQTK_WF.out.version
            FastpFull_WF( SEQTK_WF.out.reads, params.two_color_chemistry, params.adapter_sequence, params.adapter_sequence_r2,ch_publish_dir, ch_disable_publish)
            ch_fastp_report = FastpFull_WF.out.report
            ch_fastp_version = FastpFull_WF.out.version
        }
    } else {

        FastpFull_WF(ch_fastqs, params.two_color_chemistry, params.adapter_sequence, params.adapter_sequence_r2,ch_publish_dir, ch_disable_publish)
        ch_fastp_report = FastpFull_WF.out.report
        ch_fastp_version = FastpFull_WF.out.version
        
     }

    // Process fastp output and branch based on filtered read count
    FastpFull_WF.out.reads
        .join(FastpFull_WF.out.json)
        .map { sample_id, reads, json -> 
            def fastp_data = new groovy.json.JsonSlurper().parseText(json.text)
            def read_count = fastp_data.summary.after_filtering.total_reads
            [sample_id, reads, read_count]
        }
        .branch {
            small: it[2] < params.min_reads
            large: it[2] >= params.min_reads
        }
        .set { branched_reads_after_filter }
    ch_fastqs_filtered_poor_quality_samples = branched_reads_after_filter.large.map { sample_id, reads, read_count ->
            tuple(sample_id, reads)
    }
              
    SALMONQUANT (ch_fastqs_filtered_poor_quality_samples, ch_salmon_index, ch_publish_dir, ch_disable_publish)
    ch_salmon_report = SALMONQUANT.out.outdir
    ch_salmon_version = SALMONQUANT.out.version
    
    CREATE_TXIMPORT_SALMON_TX_GENE ( SALMONQUANT.out.outdir, ch_tx2gene, ch_publish_dir, ch_disable_publish)
    df_tximport = CREATE_TXIMPORT_SALMON_TX_GENE.out.collect()
    MERGE_TXIMPORT_SALMON_TX_GENE_WF (df_tximport, ch_publish_dir, ch_enable_publish)
    ch_qunt_merge_tsv = MERGE_TXIMPORT_SALMON_TX_GENE_WF.out.quant_merge_tsv
    CREATE_MATRIX_SALMON_TX_GENE ( ch_qunt_merge_tsv, ch_publish_dir, ch_enable_publish )

    STARALIGN ( ch_fastqs_filtered_poor_quality_samples, ch_star_index, ch_publish_dir,ch_disable_publish)
    ch_star_report = STARALIGN.out.outdir
    ch_star_version = STARALIGN.out.version
    ch_junction_file = STARALIGN.out.junction
    
    STARGETPRIM_WF ( STARALIGN.out.bam, ch_publish_dir, ch_enable_publish)
    ch_stargetprim_version = STARGETPRIM_WF.out.version
    ch_raw_bam = STARGETPRIM_WF.out.primary_bam.collect()

    // STARGETPRIM_WF.out.junction
    //     .collectFile( name: "junction_files.txt", newLine: true, sort: { it[0] }, storeDir: "${params.tmp_dir}" )
    //     { it[0].replaceFirst(/_.*/,"") + "\t" + "${params.publish_dir}_${params.timestamp}/secondary_analyses/alignment_htseq/" + it[1].getName() }

    // STARGETPRIM_WF.out.primary_bam
    //     .collectFile( name: "bam_files.txt", newLine: true, sort: { it[0] }, storeDir: "${params.tmp_dir}" )
    //         { it[0] + "\t" + "${params.publish_dir}_${params.timestamp}/secondary_analyses/alignment_htseq/" + it[1].getName() }

    GENE_BODY_COVERAGE_RNA (STARGETPRIM_WF.out.primary_bam, ch_genebody_ref, ch_disable_publish )
    ch_gene_body_df = GENE_BODY_COVERAGE_RNA.out.df.collect()
    GENE_BODY_COVERAGE_RNA_PLOT (ch_gene_body_df, ch_disable_publish)
    ch_gene_body_report = GENE_BODY_COVERAGE_RNA_PLOT.out.gene_body_png
               
    QUALIMAP_BAMRNA ( STARGETPRIM_WF.out.bam_raw , ch_gtf, ch_publish_dir, ch_disable_publish)
    ch_qualimap_report = QUALIMAP_BAMRNA.out.outdir
    ch_qualimap_version = QUALIMAP_BAMRNA.out.version
                        
                    
    HTSEQ_COUNTS ( STARGETPRIM_WF.out.primary_bam, ch_gtf, ch_publish_dir, ch_disable_publish)
    ch_htseq_counts_version = HTSEQ_COUNTS.out.version
    
    CREATE_HTSEQ_SUMMARY_WF ( HTSEQ_COUNTS.out.htseq_counts, ch_tx2gene, ch_publish_dir, ch_disable_publish)
    df_htseqsum = CREATE_HTSEQ_SUMMARY_WF.out.htseq_summary.collect()
    MERGE_HTSEQ_SUMMARY_WF ( df_htseqsum, ch_publish_dir, ch_enable_publish)
    ch_merge_tsv = MERGE_HTSEQ_SUMMARY_WF.out.merge_tsv.collect()
    CREATE_HTSEQ_MATRIX ( ch_merge_tsv, ch_publish_dir, ch_enable_publish)
    PLOTTER_PCAHEATMAP_HTSEQ_SUMMARY_WF (ch_merge_tsv, ch_publish_dir, ch_disable_publish)
    ch_pcaheatmap_plot=  PLOTTER_PCAHEATMAP_HTSEQ_SUMMARY_WF.out.pcaheatmap_plot
    ch_reference_celltype =  Channel.fromPath(params.celltype_ref)
    CELL_TYPING (ch_merge_tsv, ch_reference_celltype.collect(), ch_publish_dir, ch_enable_publish)
     
    //STARFUSION_WF (ch_junction_file, params.ctat_resource_lib, ch_publish_dir, ch_disable_publish)
    //ch_starfusion_all = STARFUSION_WF.out.fusion_file.collect()
    //MERGE_STAR_FUSION ( ch_starfusion_all, ch_publish_dir, ch_enable_publish )
      

    /*
    ========================================================================================
        REPORTING
    ========================================================================================
    */

    
    CREATE_QC_REPORT (ch_qualimap_report.collect(), ch_publish_dir, ch_enable_publish)
    ch_qc_report_stats= CREATE_QC_REPORT.out.stats
    ch_qc_report_pstats = CREATE_QC_REPORT.out.pstats
    
    CALC_DYNAMICRANGE (CREATE_HTSEQ_MATRIX.out.htseq_matrix ,ch_publish_dir, ch_enable_publish)

    if ( params.skip_subsampling || branched_reads.large.ifEmpty { true } ) {
    	collect_master_stats = ch_merge_tsv.collect().ifEmpty([])
	    .combine( ch_fastp_report.collect().ifEmpty([]))
	    .combine( ch_qunt_merge_tsv.collect().ifEmpty([]) )
	    .combine( ch_qc_report_stats.collect().ifEmpty([]) )
	    .combine( CELL_TYPING.out.annotations.collect().ifEmpty([]) )
	    .combine( ch_gene_body_df.collect().ifEmpty([]) )
	    .combine( CALC_DYNAMICRANGE.out.collect().ifEmpty([]) )
	    
    }else{

       collect_master_stats = ch_merge_tsv.collect().ifEmpty([])
            .combine( ch_fastp_report.collect().ifEmpty([]))
            .combine( SEQTK_WF.out.metadata.collect().ifEmpty([]))
            .combine( ch_qunt_merge_tsv.collect().ifEmpty([]) )
            .combine( ch_qc_report_stats.collect().ifEmpty([]) )
            .combine( CELL_TYPING.out.annotations.collect().ifEmpty([]) )
            .combine( ch_gene_body_df.collect().ifEmpty([]) )
            .combine( CALC_DYNAMICRANGE.out.collect().ifEmpty([]) )


    }

    CREATE_MASTER_STATS (collect_master_stats, COUNT_READS_FASTQ_WF.out.combined_read_counts, branched_reads.large.ifEmpty{}, ch_publish_dir, ch_disable_publish)

    ch_tool_versions = ch_fastp_version.take(1)
                                    .combine(ch_seqtk_version.take(1).ifEmpty([]))
                                    .combine(ch_salmon_version.take(1).ifEmpty([]))
                                    .combine(ch_star_version.take(1).ifEmpty([]))
                                    .combine(ch_qualimap_version.take(1).ifEmpty([]))
                                    .combine(ch_stargetprim_version.take(1).ifEmpty([]))
                                    .combine(ch_htseq_counts_version.take(1).ifEmpty([]))


    REPORT_VERSIONS_WF(
                            ch_tool_versions,
                            ch_publish_dir,
                            ch_enable_publish)

    collect_mqc = CREATE_MASTER_STATS.out.stats.collect().ifEmpty([])
                    .combine( ch_fastp_report.collect().ifEmpty([]))
                    .combine( ch_star_report.collect().ifEmpty([]))
                    .combine( ch_qualimap_report.collect().ifEmpty([]))
                    .combine( ch_salmon_report.collect().ifEmpty([]))
                    .combine( ch_qc_report_stats.collect().ifEmpty([]))
                    .combine( ch_qc_report_pstats.collect().ifEmpty([]))
                    .combine( ch_qunt_merge_tsv.collect().ifEmpty([]))
                    .combine( ch_merge_tsv.collect().ifEmpty([]))
                    .combine( ch_pcaheatmap_plot.collect().ifEmpty([]))
                    .combine( ch_gene_body_report.collect().ifEmpty([]))
                    .combine( REPORT_VERSIONS_WF.out.versions.collect().ifEmpty([]))

    params_meta = [
            session_id: workflow.sessionId,
            instrument: params.instrument,
            genome: params.genome,
            read_length: params.read_length,
            adapter_sequence: params.adapter_sequence,
            adapter_sequence_r2: params.adapter_sequence_r2,
            skip_subsampling: params.skip_subsampling,
            skip_10X: params.skip_10X
    ]

    if (!params.skip_subsampling) {
        params_meta['subsample_reads'] = params.n_reads
    }

    if (!params.skip_10X) {
        params_meta['lib_protocol_10x'] = params.lib_protocol_10x
        params_meta['lib_type_10x'] = params.lib_type_10x
    }

    MULTIQC_WF ( collect_mqc,
                    params_meta,
                    ch_multiqc_config.collect(),
                    ch_publish_dir,
                    ch_enable_publish
                  )
        
                            
    //ch_deliverable_files = CUSTOM_FASTQ_MERGE_WF.out.merged_fastqs_to_compress.collect().ifEmpty([])
    //                                .combine(ch_raw_bam.ifEmpty([]))
    //                                .combine(df_htseqsum.ifEmpty([]))
    //                                .combine(df_tximport.ifEmpty([]))
    //                                .combine(ch_merge_tsv.ifEmpty([]))
    //                               .combine(ch_qunt_merge_tsv.ifEmpty([]))

    //CUSTOM_PIPELINE_JSON_OUTPUT (   ch_deliverable_files,
    //                              ch_publish_dir,
    //                              ch_enable_publish
    //                   )
                                       
                                       
    
   
  
}

workflow {
    
      
    if ( params.reads ) {

        ch_reads = Channel.fromFilePairs( params.reads , size: -1 , checkExists: true )
                            .map { tag, pair -> subtags = ( tag =~ /(.*)_(S\d+)_(L+\d+)/)[0]; [ subtags[1], subtags[2], subtags[3], pair ] }
                            .groupTuple()
                            .map { tag, sample, lane, pair -> [ tag, sample.flatten(), lane.flatten(), pair.flatten() ] }

    } else if ( params.input_csv  ) {
        
        ch_reads = Channel.fromPath( params.input_csv  ).splitCsv( header:true )
                            .map { row -> [ row.sampleId, [ row.read1, row.read2 ] ] }
                            .map { tag, pair -> subtags = ( tag =~ /(.*)_(S\d+)_(L+\d+)/)[0]; [ subtags[1], subtags[2], subtags[3], pair ] }
                            .groupTuple()
                            .map { tag, sample, lane, pair -> [ tag, sample.flatten(), lane.flatten(), pair.flatten() ] }
    
    }

    ch_multiqc_config = Channel.fromPath(params.multiqc_config, checkIfExists: true)
    ch_reference_celltype =  Channel.fromPath(params.celltype_ref)
    
     RNASEQ_WF( 
                params.publish_dir,
                params.enable_publish,
                params.disable_publish,
                ch_reads, 
                params.input_csv,
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

// workflow.onComplete {
//     output                          = [:]
//     output["pipeline_run_name"]     = workflow.runName
//     output["pipeline_name"]         = workflow.manifest.name
//     output["pipeline_version"]      = workflow.manifest.version
//     output["pipeline_session_id"]   = workflow.sessionId
//     output["output"]                = [:]
//     output["output"]["bam"]         = [:]
//     // output["output"]["junction"]    = [:]


    bam_outfile = file("$params.tmp_dir/bam_files.txt")
    bam_outfile_lines = bam_outfile.exists() ? bam_outfile.readLines() : []
    for ( bam_line : bam_outfile_lines ) {
        def (sample_name, bam_path) = bam_line.split('\t')
        output["output"]["bam"][sample_name] = [:]
        output["output"]["bam"][sample_name]["bam"] = bam_path
    }
    
//     // junction_outfile = file("$params.tmp_dir/junction_files.txt")
//     // junction_outfile_lines = junction_outfile.readLines()
//     // for ( junction_line : junction_outfile_lines ) {
//     //     def (sample_name, junction_path) = junction_line.split('\t')
//     //     output["output"]["junction"][sample_name] = [:]
//     //     output["output"]["junction"][sample_name]["junction"] = junction_path
//     // }
    
//     def output_json = JsonOutput.toJson(output)
//     def output_json_pretty = JsonOutput.prettyPrint(output_json)
//     File outputfile = new File("$params.tmp_dir/output.json")
//     outputfile.write(output_json_pretty)
//     println(output_json_pretty)
// }

workflow.onError {
    output                          = [:]
    output["pipeline_run_name"]     = workflow.runName
    output["pipeline_name"]         = workflow.manifest.name
    output["pipeline_version"]      = workflow.manifest.version
    output["pipeline_session_id"]   = workflow.sessionId
    output["output"]                = [:]
    output["output"]["bam"]         = [:]
    // output["output"]["junction"]    = [:]
    
    def output_json = JsonOutput.toJson(output)
    def output_json_pretty = JsonOutput.prettyPrint(output_json)
    File outputfile = new File("$params.tmp_dir/output.json")
    outputfile.write(output_json_pretty)
    println(output_json_pretty)
    

    def subject = """\
        [nf-rnaseq-pipeline] FAILED: ${workflow.runName}
        """

    def msg = """\
        Pipeline execution summary 
        --------------------------------
        Script name       : ${workflow.scriptName ?: '-'}
        Script ID         : ${workflow.scriptId ?: '-'}
        Workflow session  : ${workflow.sessionId}
        Workflow repo     : ${workflow.repository ?: '-' }
        Workflow revision : ${workflow.repository ? "$workflow.revision ($workflow.commitId)" : '-'}
        Workflow profile  : ${workflow.profile ?: '-'}
        Workflow cmdline  : ${workflow.commandLine ?: '-'}
        Nextflow version  : ${workflow.nextflow.version}, build ${workflow.nextflow.build} (${workflow.nextflow.timestamp})
        Error Report      : ${workflow.errorReport}
        """
        .stripIndent()
    
    log.info ( msg )
    
    if ( "${params.email_on_fail}" && workflow.exitStatus != 0 ) {
        sendMail(to: "${params.email_on_fail}", subject: subject, body: msg)
    }
}
