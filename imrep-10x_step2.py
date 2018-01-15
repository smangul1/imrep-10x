import csv
import argparse
from Levenshtein import distance
import os
import sys
import collections


def pairwise_distance(l):
    d_list=[]
    for i in l:
        for j in l:
            if i>j:
                d_list.append(distance(i,j))
    return min(d_list)



ap = argparse.ArgumentParser()
ap.add_argument('inFile', help='Mapped reads in bam format')
ap.add_argument('out', help='R2 file')
args = ap.parse_args()


base=os.path.basename(args.inFile)
sample=os.path.splitext(base)[0]





#TRB,GCGGAGCTGGGATG,M02564:19:000000000-G15E1:1:2104:16904:19318,CASSSGTANQPQHF,2

file=open(args.inFile)
reader=csv.reader(file)

dict_TCR={}

barcodes=set()

for line in reader:

    
    barcodes.add(line[1])

file.close()


print ("Number barcodes", len(barcodes))

for b in barcodes:
    dict_TCR[b]=([],[])

file=open(args.inFile)
reader=csv.reader(file)
next(reader,None)

for line in reader:


    type=line[0]
    b=line[1]
    cdr3=line[3]

    if type=="TRA":
        dict_TCR[b][0].append(cdr3)
    elif type=="TRB":
        dict_TCR[b][1].append(cdr3)


file.close()

file=open(args.out+".ed","w")
file_unique=open(args.out+".unique.CDR3.csv","w")
file_multiple_CDR3=open(args.out+".multiple.CDR3.csv","w")



#file - with only barcode with >1 CDR3s
file.write("sample,chain type, barcode, CDR3s, min ed")
file.write("\n")

#file_unique
file_unique.write("sample,chain type, barcode, CDR3, nReads")
file_unique.write("\n")


#file_multiple_CDR3
file_multiple_CDR3.write("sample,chain type, barcode, CDR3, nReads")
file_multiple_CDR3.write("\n")

##sample,chain type, barcode, CDR3, number of CDR3 for this barcode

for key,value in dict_TCR.items():

    counter_TCRA = collections.Counter(value[0])
    counter_TCRB = collections.Counter(value[1])







    if len(counter_TCRA)==1:

        cdr3_TCRA=next (iter (counter_TCRA.keys()))
        count_TCRA=next (iter (counter_TCRA.values()))

        file_unique.write(sample + ",TCRA," + str(key) + ","+str(cdr3) + "," + str(count))
        file_unique.write("\n")

    if len(counter_TCRA) > 1:
        for key2, value2 in counter_TCRA.items():
            cdr3 = key2
            count = value2
            file_multiple_CDR3.write(sample + ",TCRA," + str(key) + "," + str(cdr3_TCRA) + "," + str(count_TCRA))
            file_multiple_CDR3.write("\n")

    if len(counter_TCRB) == 1:
        cdr3_TCRB = next(iter(counter_TCRB.keys()))
        count_TCRB = next(iter(counter_TCRB.values()))

        file_unique.write(sample + ",TCRB," + str(key) + "," + str(cdr3_TCRB) + "," + str(count_TCRB))
        file_unique.write("\n")

    if len(counter_TCRB) > 1:
        for key2, value2 in counter_TCRB.items():
            cdr3 = key2
            count = value2
            file_multiple_CDR3.write(sample + ",TCRB," + str(key) + "," + str(cdr3) + "," + str(count))
            file_multiple_CDR3.write("\n")

file.close()
file_unique.close()
file_multiple_CDR3.close()



	

