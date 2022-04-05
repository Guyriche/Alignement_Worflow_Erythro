nextflow.enable.dsl=2

params.reads = "data/reads/ENCSR000COQ1_1.fastq.gz"
params.ref   = "data/Refs/*.fa"
params.outdir = "results"
read_pairs_ch = Channel.fromFilePairs( params.reads )

    process INDEXING_REF {
    input :
         path genome

    output:
        path "${genome}"

    script:
        """
        samtools faidx ${genome}
        """
    }


    process ALIGN {
        publishDir "${params.outdir}/Align", mode:'copy'
        input:
            path genome
            path reads

        output:
            path "${reads}.sam"

        script:
        """
          minimap2 -ax map-ont ${genome} ${reads} > ${reads}.sam
        """
    }

workflow{
    genome_ch = channel.fromPath(params.ref, checkIfExists:true)
    reads_ch = channel.fromPath(params.reads, checkIfExists:true)

    index_ch  = INDEXING_REF(genome_ch)
    align_ch  = ALIGN(index_ch, reads_ch)
}