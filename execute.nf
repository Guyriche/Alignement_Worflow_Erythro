//nextflow.enable.dsl=2

/*
 * Le pipeline neccessite l'installation des programmes suivants:
 *  - Minimap2 (https://github.com/lh3/minimap2)
 *  - samtools (http://www.htslib.org/)
 *  - Racon
 */

/*
* Define the default parameters for input
*/

libraries = channel
    .fromPath(file('data/Refs/*.fasta'))
    .collect()
    .view()

 //libraries = [ file('data/Refs/GATA1.fasta'), file('data/Refs/GYPB_seq.fasta')  ]
 params.genome  = "data/Refs/*.fa"
 params.reads   = "data/Reads/ENCSR000COQ1_{1,2}.fastq.gz"
 params.read    = "data/Reads/FAS430170.fastq"
 params.outdir = "results"
 params.racon_nb = 4
 params.one = 1


 println """\
          ERYTHROCYTE - N F   P I P E L I N E
          ===================================
          genome       : ${params.genome}
          reads        : ${params.reads}
          read         : ${params.read}
          outdir       : ${params.outdir}
          """
          .stripIndent()


 /*
 * Process  :    Create a FastA Genome index with minimap2 or Samtools
 * input    :    genome file(FastA)
 * output   :    FastA.fai file (index)
 * Tools    :    Samtools
 */
    genome_ch = channel.fromPath(params.genome, checkIfExists:true)
    //libraries = [ file('data/Refs/genome.fa')]

    process INDEXING_REF {

    tag "Indexing the ${genome}"
    publishDir "${params.outdir}/genome_index_ch", mode:'copy'

    input :
         each path (fasta) from libraries

    output:
        path "${fasta}.fai" into genome_index_ch

    script:
         """
         samtools faidx ${fasta}
         """
    }

 /*
 * Process : align the reads with the Ref Use Minimap2 for Oxford Nanopore genomic reads  and Convert Sam -> Bam
 * Input   :
 *           Ref fastA file
 *           Reads fastq file
 * Output  : file Sam/Bam
 * Tools :   Minimap2
 */

    reads_ch = channel.fromPath(params.read, checkIfExists:true)

    process ALIGN {

        tag "Align the ${seq} in relation to ${fasta}"
        publishDir "${params.outdir}/genome_alignToSam_ch", mode:'copy'

      input:
          file seq from reads_ch
          each path(fasta) from libraries

      output:
          path "${seq.baseName}_${fasta.baseName}.sam" into genome_alignToSam_ch

      """
      minimap2 -ax map-ont ${fasta} ${seq} > "${seq.baseName}_${fasta.baseName}.sam"
      """
    }

 /*
 * Process :     Convert SAM to BAM
 * Input   :     File SAM
 * Output  :     File Bam which is an indexed file
 * Tools   :     Samtools
 */

    process CONVERT_SAM_TO_BAM {

        tag "CONVERT ${file_sam} TO ${file_sam}.bam"
        publishDir "${params.outdir}/filesBam_ch", mode:'copy'

        input:
            path file_sam from genome_alignToSam_ch

        output:
             path "${file_sam.baseName}.bam" into filesBam_ch, filesBam2_ch, filesBam3_ch, filesBam4_ch

        script:
            """
            samtools view -bh ${file_sam} > ${file_sam.baseName}.bam
            samtools sort ${file_sam.baseName}.bam -o ${file_sam.baseName}.bam
            """
    }


/*
* Process :     Indexing of the bam file in order to be able to perform a visualization on IGV
* Input   :     File Bam
* Output  :     File Bam.bai which is an indexed file
* Tools   :     Samtools
*/

  process INDEXING_BAM {

        tag "Indexing ${file_bam}"
        publishDir "${params.outdir}/filesBam_ch", mode:'copy'

        input:
            path file_bam from filesBam_ch

        output:
             path "${file_bam}.bai" into bam_Index_ch

        script:
            """
            samtools index -b ${file_bam} > ${file_bam}.bai
            """
  }

 /*
  * Process :   Polissage With RACON to remove errors that can persist through to the novo assemblies
  * Input   :
  *             File Sam
  *             File Fasta indexed
  *             File fastq (reads)
  * Output  :   File Fasta
  * Tools   :   Racon
  */

    readsForPolissage_ch = channel.fromPath(params.read, checkIfExists:true)

    process POLISSAGE {

        tag "quantification on ${fasta}"
        publishDir "${params.outdir}/racon_fasta_ch", mode:'copy'

        input:
            each path(fasta) from libraries
            path seqs from readsForPolissage_ch

        output:
            path "${seqs.baseName}_${fasta.baseName}${params.racon_nb}.fasta" into racon_fasta_ch, racon_consesnsus

        script:
        """
        set  +eu
        subp=all
        for i in `seq 1 ${params.racon_nb}`; do
            ii=\$((\$i-1))
            if [ "\$i" == "${params.one}" ]
            then
                samtools faidx ${fasta}
                minimap2 -ax map-ont ${fasta} ${seqs} > ${seqs.baseName}_${fasta.baseName}\$i.sam
                racon -m 8 -x -6 -g -8 -w 500 -t 7 ${seqs} ${seqs.baseName}_${fasta.baseName}\$i.sam ${fasta} > ${seqs.baseName}_${fasta.baseName}\$i.fasta
            else
                samtools faidx ${seqs.baseName}_${fasta.baseName}\$ii.fasta
                minimap2 -ax map-ont ${seqs.baseName}_${fasta.baseName}\$ii.fasta ${seqs} > ${seqs.baseName}_${fasta.baseName}\$i.sam
                racon -m 8 -x -6 -g -8 -w 500 -t 7 ${seqs} ${seqs.baseName}_${fasta.baseName}\$i.sam ${seqs.baseName}_${fasta.baseName}\$ii.fasta > ${seqs.baseName}_${fasta.baseName}\$i.fasta
            fi
        done

        """
    }

/*
* Process :     Identification of Reads variants on the REF Genome
* Input   :
*               Bam File returned by Racon indexed
*               Fasta file returned by Racon
* Output  :
*               VCF File
* Tools :       BcfTools
*/


    process PRINTVCF {

        tag "Polishing on ${genome}"
        publishDir "${params.outdir}/vcfFile", mode:'copy'

        input:
            path genome from racon_fasta_ch
            path files_bam from filesBam2_ch

        output:
            path "${genome.baseName}.vcf" into vcfFile, vcfFile2, vcfFile3

        script:
        """
        bcftools mpileup -a 'INFO/AD' -Ov -f ${genome} ${files_bam} | bcftools call -mv -o ${genome.baseName}.vcf
        """
    }


/*
* Process :     Modified VCF file to contains ALT(%) and REF(%)
* Input   :
*               Bam File returned by Racon indexed
*               VCF file
* Output  :
*               VCF File(Add Extension CHU)
* Tools :       Script Python
*/

process MODIFEDVCF {

    tag "Modified vcf ${fileBam}"
    publishDir "${params.outdir}/vcfFile", mode:'copy'

    input:
        path fileBam from filesBam3_ch
        path fileVcf from vcfFile



    script:
    """
    python3 ${projectDir}/scripts/modifiedVCF.py ${fileBam} ${fileVcf} ${projectDir}/
    """
}

/*
* Process :     Comparison between the ref Genome and the consensus file
* Input   :
*               Fasta File REF
*               Fasta file returned by Racon
* Output  :
*               Delta File
* Tools :       Mummer (nucmer)
*/

/*

process MUMMER {

    tag "Nucmer on ${fasta}"
    publishDir "${params.outdir}/nucmerFile", mode:'copy'
    input:
        each path (fasta) from libraries
        path consensus from racon_consesnsus

    output:
        path "${consensus}.snps" into nucmerFile

    script:
    """
    nucmer --prefix=SNPsChr ${fasta} ${consensus}
    show-snps -Clr SNPsChr.delta > ${consensus}.snps
    """
}

 */


/*
* Process :     Generate a consensus file
* Input   :
*               Bam File returned by Racon indexed
*               Ref fasta file
* Output  :
*               Fasta File()
* Tools :       Script Python
*/

process CONSENSUS {

    tag "Generate consensus ${fileBam}"
    publishDir "${params.outdir}/consensus", mode:'copy'

    input:
        path fileBam from filesBam4_ch
        each path (fasta) from libraries

    script:
    """
    python3 ${projectDir}/scripts/generateConsensus.py ${fileBam} ${fasta}
    """
}

process GET_STATS_VCF {

    tag "Generate consensus ${vcffile}"
    publishDir "${params.outdir}/statsVcf", mode:'copy'

    input:
        path vcffile from vcfFile2

    output:
        path "stats.txt" into statsVcf

    script:
    """
    vcf-stats ${vcffile} > stats.txt
    """

}