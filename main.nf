nextflow.enable.dsl=2

/*
 * Le pipeline neccessite l'installation des programmes suivants:
 *  - Minimap2 (https://github.com/lh3/minimap2)
 *  - samtools (http://www.htslib.org/)
 *  - Racon
 */

/*
* Define the default parameters for input
*/

 params.genome  = "data/genome.fa"
 params.reads   = "data/reads/ENCSR000COQ1_{1,2}.fastq.gz"
 params.read    = "data/reads/ENCSR000COQ1_1.fastq.gz"
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

/*
* Process : align the reads with the Ref Use Minimap2 for Oxford Nanopore genomic reads  and Convert Sam -> Bam
* Input   :
*           Ref fastA file
*           Reads fastq file
* Output  : file Sam/Bam
* Tools :   Minimap2
*/

    process ALIGN {
        tag "Align the ${read_id}"
        publishDir "${params.outdir}/Align", mode:'copy'
    	input:
    	    path genome
    	    tuple val(read_id), path(reads)

    	output:
    	    path "${read_id}.sam"

    	script:
        """
          minimap2 -ax map-ont ${genome} ${reads[0]} > ${read_id}.sam
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
        publishDir "${params.outdir}/FilesBam", mode:'copy'
        input:
            path file_sam

        output:
             path "${file_sam}.bam"

        script:
            """
            samtools view -bh ${file_sam} > ${file_sam}.bam
            samtools sort ${file_sam}.bam -o ${file_sam}.bam
            """
  }

/*
* Process :     Indexing of the bam file in order to be able to perform a visualization on IGV
* Input   :     File Bam
* Output  :     File Bam.bai which is an indexed file
* Tools   :     Samtools
*/

  process INDEXING_BAM {
        tag "Indexing ${file_sam}"
        publishDir "${params.outdir}/FilesBam", mode:'copy'
        input:
            path file_sam

        output:
             path "${file_sam}.bai"

        script:
            """
            samtools index -b ${file_sam} > ${file_sam}.bai
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

    process POLISSAGE {
        tag "quantification on ${sample_id}"
        publishDir "${params.outdir}/Racon", mode:'copy'
        input:
            path genome
            tuple val(sample_id), path(reads)

        output:
            path "${sample_id}${params.racon_nb}.fasta"
            //path "${sample_id}${params.racon_nb}.sam"

        script:
        """
        set  +eu
        ln -s ${genome} ${sample_id}0.fasta
        for i in `seq 1 ${params.racon_nb}`; do
            ii=\$((\$i-1))
            samtools faidx ${sample_id}\$ii.fasta
            minimap2 -ax map-ont ${sample_id}\$ii.fasta ${reads[0]} > ${sample_id}\$i.sam
            racon -m 8 -x -6 -g -8 -w 500 -t 7 ${reads[0]} ${sample_id}\$i.sam ${sample_id}\$ii.fasta > ${sample_id}\$i.fasta
        done
        """
    }

/*
* Process : Script python to compare two File
* Input   : 02 File to compare
* Output  : Boolean value (True if FileA eq FileB ? False)
*/

  process COMAPREFILE {

        input:
            path file1
            path file2

        script:
        """
        #!/usr/bin/env python3

        import filecmp

        value = filecmp.cmp('${file1}', '${file2}')
        print(value)
        """
  }

  process ITERATION_OF_POLISSAGE {

    tag "Polishing on ${sample_id}"
    publishDir "${params.outdir}/Racon", mode:'copy'
    input:
        path genome
        tuple val(sample_id), path(reads)

    output:
        path "${sample_id}_${params.racon_nb}.fasta"
        path "${sample_id}_${params.racon_nb}.sam"

    script:
    """
    set  +eu
    subp=all
    for i in `seq 1 ${params.racon_nb}`; do
        ii=\$((\$i-1))
        if [ "\$i" == "${params.one}" ]
        then
            samtools faidx ${genome}
            minimap2 -ax map-ont ${genome} ${reads[0]} > ${sample_id}\$i.sam
            racon -m 8 -x -6 -g -8 -w 500 -t 7 ${reads[0]} ${sample_id}\$i.sam ${genome} > ${sample_id}\$i.fasta
        else
            samtools faidx ${sample_id}\$ii.fasta
            minimap2 -ax map-ont ${sample_id}\$ii.fasta ${reads[0]} > ${sample_id}\$i.sam
            racon -m 8 -x -6 -g -8 -w 500 -t 7 ${reads[0]} ${sample_id}\$i.sam ${sample_id}\$ii.fasta > ${sample_id}\$i.fasta
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
    publishDir "${params.outdir}/VcfFile", mode:'copy'

    input:
        path genome
        path file_bam

    output:
        path "${genome}.vcf"

    script:
    """
    bcftools mpileup -Ov -f ${genome} ${file_bam} | bcftools call -mv -o ${genome}.vcf
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

process MUMMER {

    tag "Nucmer on ${genome}"
    publishDir "${params.outdir}/NucmerFile", mode:'copy'
    input:
        path genome
        path consensus

    output:
        path "${consensus}.snps"

    script:
    """
    nucmer --prefix=SNPsChr ${genome} ${consensus}
    show-snps -Clr SNPsChr.delta > ${consensus}.snps
    """
}

/*
* Process :     Identifying Heterozygous and Homozygous Variants in VCF File
* Input   :
*               VCF File
* Output  :
*               VCF File
* Tools :       VcfTools
*/

process MUTATION {

    tag "Mutation on ${vcfFile}"
    publishDir "${params.outdir}/Mutation", mode:'copy'

    input:
        path vcfFile
    output:
        path "${vcfFile}_out"

    script:
    """
    vcftools --vcf ${vcfFile} --freq --out ${vcfFile}_out
    """
}

process MODIFEDVCF {

    tag "Modified vcf ${fileBam}"
    publishDir "${params.outdir}/VcfFile", mode:'copy'

    input:
        path fileBam
        path newVcf

    output:
        path "${newVcf}"

    script:
    """
    python3 ${projectDir}/scripts/modifiedVCF.py ${fileBam} ${newVcf}
    """
}

process CONSENSUS {

    tag "Generate consensus ${fileBam}"
    publishDir "${params.outdir}/consensus", mode:'copy'

    input:
        path fileBam
        path genome

    script:
    """
    python3 ${projectDir}/scripts/generateConsensus.py ${fileBam} ${genome}
    """
}

process COPYFILE {

    publishDir "${params.outdir}/VcfFile", mode:'copy'

    input:
        path newVcf

    output:
        path "${newVcf}"
    script:
    """
    cp -r ${newVcf} ${projectDir}/results/VcfFile/
    """

}

workflow REINDEXING_BAM{
    take: data
    main:
        INDEXING_BAM(data)
    emit:
        INDEXING_BAM.out
}

workflow REALIGN_FASTA{
    take:
        fasta
        reads
    main:
        ALIGN(fasta, reads)
    emit:
        ALIGN.out
}

workflow RECONVERT{
    take: data
    main:
        CONVERT_SAM_TO_BAM(data)
    emit:
        CONVERT_SAM_TO_BAM.out
}

workflow {

    genome_ch = channel.fromPath(params.genome, checkIfExists:true)
    reads_ch = channel.fromFilePairs(params.reads, checkIfExists:true)
    read_ch     = channel.fromPath(params.read, checkIfExists:true)

    index_ch    = INDEXING_REF(genome_ch)
    align_ch    = ALIGN(index_ch, reads_ch)
    convert_sam_to_bam_ch = CONVERT_SAM_TO_BAM(align_ch)
    indexBam_ch = INDEXING_BAM(convert_sam_to_bam_ch)
    racon_ch    = POLISSAGE(genome_ch, reads_ch)
    vcfFile_ch = PRINTVCF(racon_ch, convert_sam_to_bam_ch)
    newVcf_ch  = MODIFEDVCF(convert_sam_to_bam_ch, vcfFile_ch)
    COPYFILE(newVcf_ch)
    nucmer_ch   = MUMMER(genome_ch, racon_ch)
    consensus_ch= CONSENSUS(convert_sam_to_bam_ch, genome_ch)
    //mutation_ch = MUTATION(vcfFile_ch)

}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/align_report.html\n" : "Oops .. something went wrong" )
}

// nextflow run main.nf -resume -with-report -with-trace -with-timeline -with-dag dag.png
// cmp $1 $2 1>/dev/null 2>&1; resultat=$?
