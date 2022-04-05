import pysam
import matplotlib.pyplot as plt
import sys

input=sys.argv[1]
input2=sys.argv[2]
out = sys.argv[3]

from fuc import pyvcf
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from vcf_parser import VCFParser


def countpourcentbases(base_count_dict):
    bases = []
    pourcentBases = {}
    sumBasesOccurences = 0

    for base in ['A', 'C', 'G', 'T']:
        assert base in base_count_dict, "Wrong dictionnary in argument"
        if base_count_dict[base] > 0:
            sumBasesOccurences += base_count_dict[base]
            bases.append(base)

    for base in bases:
        pourcentBases[base] = str(int(((base_count_dict[base] / sumBasesOccurences) * 100) // 1))+"%"

    return pourcentBases


def isHomoHetero(base_count_dict):
    """
    Return true if the % of the base that are the most accurate is between 40 et 70%
    :param base_count_dict: {'A': x, 'G': y, etc} where keys are bases and x, y, etc the occurence of the base
    :return:
    """
    bases = ['A', 'C', 'G', 'T']
    # Number of bases where value is > 0
    nbBases = 0
    # Sum of bases occurences
    sumBasesOccurences = 0

    for base in bases:
        assert base in base_count_dict, "Wrong dictionnary in argument"
        # Check that the value of base in dict is higher than 0
        if base_count_dict[base] > 0:
            nbBases += 1
            sumBasesOccurences += base_count_dict[base]

    if nbBases >= 2:
        # Get key that has the maximum value
        find_max = max(base_count_dict, key=base_count_dict.get)
        # Compute its % based on the sumBasesOccurences
        pourcentage_max = base_count_dict[find_max] / sumBasesOccurences
        if 0.4 < pourcentage_max < 0.7:
            return True

    return False


def processBamFile(path):
    """
    Process the Bam file
    :param path: the path to the bam file
    :return: return a dict where each KEY is the position of the seq where all reads represents homozygote or
    heterozygote. VALUE is a dictionnary where each bases that occurs in the stack of reads has as value its coverage
    in %
    """
    # Read file
    print(path)
    samfile = pysam.AlignmentFile(path, "rb")
    total_size = 0
    # Use for the plot
    positions = []
    coverage = []
    # Store the return value
    vcf_id = {}

    # One pileupcolumn for each base in ref genome (.fasta)
    for pileupcolumn in samfile.pileup():
        total_size += 1
        # Build the lists for the plot
        positions.append(pileupcolumn.reference_pos)
        coverage.append(pileupcolumn.nsegments)

        # print("\ncoverage at base %s (index of base in ref) = %s (number of reads mapping this column)" %
        #     (pileupcolumn.reference_pos, pileupcolumn.nsegments))

        # Count the bases for a column
        base_count = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
        # All read align to this column
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # Increment the bases in base_count by 1
                if pileupread.alignment.query_sequence[pileupread.query_position] != 'N':
                    base_count[pileupread.alignment.query_sequence[pileupread.query_position]] += 1

        # Check repartition of the bases for this column
        if isHomoHetero(base_count):
#             print("\ncoverage at base %s (index of base in ref) = %s (number of reads mapping this column)" %
#                   (pileupcolumn.reference_pos, pileupcolumn.nsegments))
#             print(base_count)
#             print(countpourcentbases(base_count))
#             print(pileupcolumn.reference_name)
            vcf_id[pileupcolumn.reference_pos + 1] = countpourcentbases(base_count)

        # Change value of special_pos to get info for special position
        # # CODE Start here
        # special_pos = 548
        # if pileupcolumn.reference_pos == special_pos - 1:
        #     print(f"\n# Special position {special_pos} is {base_count}")
        # # CODE End here

    print(f"Total pileupcolumns is {total_size}")
    # plt.plot(positions, coverage)
    # plt.show()
    samfile.close()
    return vcf_id


def vcf_id_convertor(dic):
    s = ''
    for key in dic:
        s += f'{key}{dic[key]};'
    return s


def modifyVCF(path, vcf_id):
    print("Writing VCF file..")
    vf = pyvcf.VcfFrame.from_file(path)
    data = vf.df.values

    for position in vcf_id:
        for i in range(len(data)):
            if data[i][1] == position:
                vf.df.at[i, 'ALT'] = vcf_id_convertor(vcf_id[position])

    if path[-4:] == '.vcf':
        path = f'{path[:-4]}.CHU.vcf'
    else:
        path = f'{path}.CHU.vcf'
    vf.to_file(path)


# Bam and vcf file
vcf_id = processBamFile(input)
modifyVCF(input2, vcf_id)

# mainBam()


