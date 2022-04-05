import pysam
import matplotlib.pyplot as plt
from fuc import pyvcf
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from vcf_parser import VCFParser

import sys

input=sys.argv[1]
input1 = sys.argv[2]

# New fasta file for assembly
def assemblyReads(bampath, fastapath):
    print(bampath)
    print(fastapath)
    file_in = fastapath
    file_out = f'{fastapath[:-6]}.CHU.fasta'

    with open(file_out, 'w') as f_out:
        for seq_record in SeqIO.parse(open(file_in, mode='r'), 'fasta'):
            # remove .id from .description record (remove all before first space)
            seq_record.description = ' '.join(seq_record.description.split()[1:])
            # do something (print or edit seq_record)
            print('SequenceID = ' + seq_record.id)
            print('Description = ' + seq_record.description + '\n')
            print("seq is " + seq_record.seq)
            new_seq = ""
            # Read file
            samfile = pysam.AlignmentFile(bampath, "rb")

            pileupcolumnlist = []
            pileupreadlist = []
            for pileupcolumn in samfile.pileup():
                pileupcolumnlist.append(pileupcolumn)
                pileupreadlist.append(pileupcolumn.pileups)



            start = 0
            weAt = 1
            for position in range(len(seq_record.seq)):
                # if position == 150:
                #     break
                # print(f"{position + 1} {new_seq}")

                # One pileupcolumn for each base in ref genome (.fasta)
                if start == len(pileupcolumnlist):
                    new_seq += "-"

                for i in range(start, len(pileupcolumnlist)):
                    pileupcolumn = pileupcolumnlist[i]

                    # print(f"{position} {pileupcolumn.reference_pos}")
                    if pileupcolumn.reference_pos > position+1:
                        if weAt == position+1:
                            weAt += 1
                            new_seq += "-"
                        break

                    elif pileupcolumn.reference_pos == position+1:
                        weAt += 1
                        start += 1
                        base_count = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
                        # All read align to this column

                        for pileupread in pileupreadlist[i]:
                            if not pileupread.is_del and not pileupread.is_refskip:
                                # Increment the bases in base_count by 1
                                if pileupread.alignment.query_sequence[pileupread.query_position] != 'N':
                                    base_count[pileupread.alignment.query_sequence[pileupread.query_position]] += 1

                        new_seq += max(base_count, key=base_count.get)
                        #print(f"{position} {base_count} {max(base_count, key=base_count.get)}")



            seq_record.seq = Seq(new_seq)
            # write new fasta file
            r = SeqIO.write(seq_record, f_out, 'fasta')
            if r != 1:
                print('Error while writing sequence:  ' + seq_record.id)


assemblyReads(input, input1)

