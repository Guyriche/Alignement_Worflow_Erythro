import sys

input=sys.argv[1]
input1=sys.argv[2]

def processBamFile(a):
    print(a)
    return a
def modifyVCF(b,c):
    return c+b

vcf_id = processBamFile(input)
modifyVCF(input1 ,vcf_id)
print(modifyVCF(input1 ,vcf_id))