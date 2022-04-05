params.genome  = "data/*.fa"
params.ref   = "data/Refs/*.fa"
params.read    = "data/reads/ENCSR000COQ1_1.fastq.gz"
params.outdir = "results"

sequences = Channel.fromPath('*.fa')
reads_ch = channel.fromPath(params.read, checkIfExists:true)
methods = ['regular', 'expresso']
libraries = [ file('data/Refs/genome.fa'), file('data/Refs/genome1.fa') ]


process INDEXING_REF {
    publishDir "${params.outdir}/Index_fasta", mode:'copy'
    input :
     each file(fasta) from libraries

    output:
    path "${fasta}"

    script:
    """
    samtools faidx ${fasta}
    """
}

process ALIGN {
    publishDir "${params.outdir}/toto", mode:'copy'
  input:
  file seq from reads_ch
  each file(fasta) from libraries

  output:
      path "${seq}_${fasta}.sam" into toto

  """
  minimap2 -ax map-ont ${fasta} ${seq} > "${seq}_${fasta}.sam"
  """
}

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