import csv
import argparse
from Levenshtein import distance
import os

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


print "Number barcodes", len(barcodes)

for b in barcodes:
    dict_TCR[b]=([],[])

file=open(args.inFile)
reader=csv.reader(file)
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
file2=open(args.out,"w")
file3=open(args.out+".2.csv","w")


#file - with only barcode with >1 CDR3s
file.write("sample,chain type, barcode, CDR3s, min ed")
file.write("\n")

#file2 - CDR3 are presented in filed separated by ;
file2.write("sample,chain type, barcode, CDR3s, number of CDR3 for this barcode")
file2.write("\n")


#file3 - one CDR3 per line
file3.write("sample,chain type, barcode, CDR3, number of CDR3 for this barcode")
file3.write("\n")

##sample,chain type, barcode, CDR3, number of CDR3 for this barcode

for key,value in dict_TCR.iteritems():
    if len(set(value[0]))>=1:
        list_temp=list(set(value[0]))
        tmp=";".join(list(set(value[0])))
        
        if len(set(value[0]))>1:
            file.write(sample+",TCRA"+","+str(key)+","+tmp +"," + str(pairwise_distance(list_temp)) )
            file.write("\n")
    
    
        file2.write(sample+",TCRA,"+str(key)+","+tmp+","+str(len(list_temp)))
        file2.write("\n")

        for i in list_temp:
            file3.write(sample+",TCRA,"+str(key)+","+str(i)+","+str(len(list_temp)))
            file3.write("\n")

    
    if len(set(value[1]))>=1:
        list_temp=list(set(value[1]))
        tmp=";".join(list(set(value[1])))
        
        if len(set(value[1]))>1:
            file.write(sample+",TCRB"+","+str(key)+","+tmp +"," + str(pairwise_distance(list_temp)) )
            file.write("\n")

        file2.write(sample+",TCRB,"+str(key)+","+tmp+","+str(len(list_temp)))
        file2.write("\n")

        for i in list_temp:
            file3.write(sample+",TCRB,"+str(key)+","+str(i)+","+str(len(list_temp)))
            file3.write("\n")





file.close()
file2.close()
file2.close()

	

