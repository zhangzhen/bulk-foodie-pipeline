/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                                 } from '../modules/nf-core/fastqc/main'
include { TRIMGALORE                             } from '../modules/nf-core/trimgalore/main'
include { BISMARK_ALIGN                          } from '../modules/nf-core/bismark/align/main'
include { SAMTOOLS_SORT                          } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX                         } from '../modules/nf-core/samtools/index/main'
include { BISMARK_DEDUPLICATE                    } from '../modules/nf-core/bismark/deduplicate/main'
include { METHYLEXTRACT                          } from '../modules/local/methylextract/main'
include { CALLTRACK                              } from '../modules/local/calltrack/main'
include { BEDTOOLS_BAMTOBED                      } from '../modules/nf-core/bedtools/bamtobed/main'
include { POSTPROCESSBED                         } from '../modules/local/postprocessbed/main'
include { BEDTOOLS_GENOMECOV                     } from '../modules/nf-core/bedtools/genomecov/main'
include { UCSC_BEDGRAPHTOBIGWIG                  } from '../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { MACS2_CALLPEAK                         } from '../modules/nf-core/macs2/callpeak/main'
include { CALLRATIOANDDEPTH                      } from '../modules/local/callratioanddepth/main'
include { IGVTOOLS_TOTDF as IGVTOOLS_TOTDF_RATIO } from '../modules/local/igvtools/totdf/main'
include { IGVTOOLS_TOTDF as IGVTOOLS_TOTDF_DEPTH } from '../modules/local/igvtools/totdf/main'
include { HIGHSCOREPEAKS                         } from '../modules/local/highscorepeaks/main'
include { PLOT                                   } from '../modules/local/plot/main'
include { BEDTOOLS_INTERSECT                     } from '../modules/nf-core/bedtools/intersect/main'
include { PLOTTSS                                } from '../modules/local/plottss/main'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMultiqc                   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                 } from '../subworkflows/local/utils_nfcore_bulkfoodiepipeline_pipeline'
include { FOOTPRINTING                           } from '../subworkflows/local/footprinting'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BULKFOODIEPIPELINE {
    take:
    ch_samplesheet   // channel: samplesheet read in from --input
    ch_fasta         // channel: fasta file
    ch_sizes         // channel: sizes file
    ch_bismark_index // channel: [ path(bismark index)   ]
    genome_id
    macs2_gsize
    tss
    depth
    scripts_dir

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // SEQKIT_SPLIT2(ch_samplesheet)>
    // ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    // ch_fastq_pieces = SEQKIT_SPLIT2.out.reads.flatMap { meta, chunk_list ->
    //     // Group chunks by their part number
    //     def grouped_chunks = [:]
    //     chunk_list.each { chunk ->
    //         // Extract part number from filename (handle both fastq.gz and fq.gz)
    //         def part_match = chunk.name =~ /.*_(\d+)\.(fastq|fq)\.gz$/
    //         if (part_match) {
    //             def part_num = part_match[0][1]
    //             if (!grouped_chunks.containsKey(part_num)) {
    //                 grouped_chunks[part_num] = []
    //             }
    //             grouped_chunks[part_num] << chunk
    //         }
    //     }

    //     // Return paired chunks as tuples
    //     return grouped_chunks.collect { _part_num, chunks ->
    //         tuple(meta, chunks)
    //     }
    // }

    FASTQC(ch_samplesheet)
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    TRIMGALORE(ch_samplesheet)
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions)

    BISMARK_ALIGN(
        TRIMGALORE.out.reads,
        [[:], file(ch_fasta, checkIfExists: true)],
        [[:], file(ch_bismark_index, checkIfExists: true)],
    )
    ch_alignments = BISMARK_ALIGN.out.bam
    ch_alignment_reports = BISMARK_ALIGN.out.report.map { meta, report -> [meta, report, []] }
    ch_versions = ch_versions.mix(BISMARK_ALIGN.out.versions)

    /*
    * Run deduplicate_bismark
    */
    BISMARK_DEDUPLICATE(
        BISMARK_ALIGN.out.bam
    )
    ch_alignments = BISMARK_DEDUPLICATE.out.bam
    ch_alignment_reports = BISMARK_ALIGN.out.report.join(BISMARK_DEDUPLICATE.out.report)
    ch_versions = ch_versions.mix(BISMARK_DEDUPLICATE.out.versions)

    /*
     * MODULE: Run samtools sort on aligned or deduplicated bam
     */
    SAMTOOLS_SORT(
        ch_alignments,
        [[:], []],
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    /*
     * MODULE: Run samtools index on aligned or deduplicated bam
     */
    SAMTOOLS_INDEX(
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    ch_bam_sorted_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)

    METHYLEXTRACT(
        ch_bam_sorted_bai
    )
    ch_versions = ch_versions.mix(METHYLEXTRACT.out.versions)

    CALLTRACK(
        METHYLEXTRACT.out.bed
    )
    ch_versions = ch_versions.mix(CALLTRACK.out.versions)

    BEDTOOLS_BAMTOBED(
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

    POSTPROCESSBED(
        BEDTOOLS_BAMTOBED.out.bed
    )

    BEDTOOLS_GENOMECOV(
        POSTPROCESSBED.out.frag_bed.combine(Channel.from(1)),
        ch_sizes,
        'bedGraph',
        false,
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)

    UCSC_BEDGRAPHTOBIGWIG(
        BEDTOOLS_GENOMECOV.out.genomecov,
        ch_sizes,
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions)

    ch_bedpe = POSTPROCESSBED.out.frag_bed.map { meta, bedpe ->
        [meta, bedpe, []]
    }
    MACS2_CALLPEAK(
        ch_bedpe,
        macs2_gsize,
    )
    ch_versions = ch_versions.mix(MACS2_CALLPEAK.out.versions)

    HIGHSCOREPEAKS(
        MACS2_CALLPEAK.out.peak
    )

    PLOT(
        METHYLEXTRACT.out.txt
    )

    CALLRATIOANDDEPTH(
        CALLTRACK.out.summary_score_bed
    )

    IGVTOOLS_TOTDF_DEPTH(
        CALLRATIOANDDEPTH.out.depth,
        genome_id,
        'igv',
    )

    IGVTOOLS_TOTDF_RATIO(
        CALLRATIOANDDEPTH.out.ratio,
        genome_id,
        'igv',
    )

    BEDTOOLS_INTERSECT(
        CALLRATIOANDDEPTH.out.sites.map { meta, sites ->
            [meta, sites, file(tss, checkIfExists: true)]
        },
        [[:], ch_sizes],
    )

    FOOTPRINTING(CALLRATIOANDDEPTH.out.sites, HIGHSCOREPEAKS.out.bed, depth, scripts_dir)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'bulkfoodiepipeline_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
