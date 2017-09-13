import pysam
import argparse
import csv
import sys
import os
from Bio import SeqIO



ap = argparse.ArgumentParser()
ap.add_argument('barcodes', help='Mapped reads in bam format')
ap.add_argument('f2', help='Mapped reads in bam format')
ap.add_argument('imrep_file', help='Mapped reads in bam format')
ap.add_argument('out', help='Mapped reads in bam format')


args = ap.parse_args()

base=os.path.basename(args.f2)
base2=os.path.splitext(base)[0]

f = open(args.barcodes)
csv_f = csv.reader(f)


bSet=set()
bSet2=set()

for row in csv_f:
    b=row[0].split("-")[0]
    bSet.add(b)
    bSet2.add(b[1:])
    




print "Number of barcodes", len(bSet), len(bSet2)
k_full=0
k_full_match=0
k_N=0
k_N_match=0

barcodesSet=set()

#CDR3
#M02564:19:000000000-G15E1:1:2102:15042:11150    CILGSNNNDMRF    TRA     TRAV26  NA      TRAV26-2*01:20:2        TRAJ43*01:17:0  1       1       1


dict_reads_CDR3_TRA={}
dict_reads_CDR3_TRB={}

reads_TRA=set()
reads_TRB=set()

f=open(args.imrep_file,"r")
csv_f=csv.reader(f,delimiter="\t")
next(csv_f, None)  # skip the headers
    

for row in csv_f:
    read=row[0]
    cdr3=row[1]
    if row[2]=="TRA":
        dict_reads_CDR3_TRA[read]=cdr3
        reads_TRA.add(read)
    elif row[2]=="TRB":
         dict_reads_CDR3_TRB[read]=cdr3
         reads_TRB.add(read)

    
print "Number of TRA reads", len(reads_TRA)
print "Number of TRB reads", len(reads_TRB)


print "Group CDR3s by barcodes ..."

dict_barcode_TRA={}
set_barcode_TRA=set()
dict_barcode_TRB={}
set_barcode_TRB=set()

for b in bSet2:
    dict_barcode_TRA[b]=[]
    dict_barcode_TRB[b]=[]

with open(args.f2, "rU") as handle:
    for record in SeqIO.parse(handle, "fastq") :
        barcode=record.seq
        barcodeShort=record.seq[1:]



        read=record.id
        
        if "N" in record.seq:
            if record.seq[0]=="N":
                barcode=record.seq[1:]
                if barcode in bSet2:
                    barcodesSet.add(barcode)
                    k_N_match+=1
                    if read in reads_TRA:
                        cdr3=dict_reads_CDR3_TRA[read]
                        dict_barcode_TRA[barcodeShort].append(cdr3)
                        set_barcode_TRA.add(barcodeShort)
                            #print barcodeShort,cdr3,read
                            #sys.exit(1)
                    elif read in reads_TRB:
                        cdr3=dict_reads_CDR3_TRB[read]
                        dict_barcode_TRB[barcodeShort].append(cdr3)
                        set_barcode_TRB.add(barcodeShort)




                k_N+=1
        
            else:
                print "WARNING"
                print record.id, record.seq
                print "EXIT!"
                sys.exit(1)

        else:
            k_full+=1
            if barcode in bSet:
                barcodesSet.add(barcode[1:])
                k_full_match+=1
                if read in reads_TRA:
                    cdr3=dict_reads_CDR3_TRA[read]
                    dict_barcode_TRA[barcodeShort].append(cdr3)
                    set_barcode_TRA.add(barcodeShort)

                elif read in reads_TRB:
                    cdr3=dict_reads_CDR3_TRB[read]
                    dict_barcode_TRB[barcodeShort].append(cdr3)
                    set_barcode_TRB.add(barcodeShort)







print "Number of reads with full barcode",k_full
print "Number of reads with full barcode mathing",k_full_match
print "Number of reads with not-full barcode",k_N
print "Number of reads with not-full barcode mathing",k_N_match



nCDR3_TRA=[]
nCDR3_TRB=[]

file=open(args.out+"_raw.txt","w")

for key, value in dict_barcode_TRA.iteritems():
    if len(value)>0:
        for s in set(value):
            file.write(base2 +","+"TCRA,"+key +","+s+","+str(len(set(value))) )
            file.write("\n")
        nCDR3_TRA.append(len(set(value)))


for key, value in dict_barcode_TRB.iteritems():
    if len(value)>0:
        for s in set(value):
            file.write(base2 +","+"TCRB,"+key +","+s+","+str(len(set(value))) )
            file.write("\n")
        nCDR3_TRB.append(len(set(value)))

file.close()

#TRA
file=open(args.out+"_TCR_summary.txt","w")
file.write("Sample, total_barcodes, intersection, union,  Total_TRA, TRA_1, TRA_2, TRA>=2,Total_TRB, TRB_1, TRB_2, TRB>=2, Number of reads with full barcode, Number of reads with full barcode mathing, Number of reads with not-full barcode, Number of reads with not-full barcode mathing")
file.write("\n")
total_barcodes=len(barcodesSet)
Total_TRA=len(nCDR3_TRA)
TRA_1=nCDR3_TRA.count(1)
TRA_2=nCDR3_TRA.count(2)
Total_TRB=len(nCDR3_TRB)
TRB_1=nCDR3_TRB.count(1)
TRB_2=nCDR3_TRB.count(2)
TRB_gt2=len(nCDR3_TRB) - nCDR3_TRB.count(1) -nCDR3_TRB.count(2)




intersection=len(set_barcode_TRA.intersection(set_barcode_TRB))
union=len( set_barcode_TRA | set_barcode_TRB)




TRA_gt2=len(nCDR3_TRA) - nCDR3_TRA.count(1) -nCDR3_TRA.count(2)



file.write(base2 +"," + str(total_barcodes) +  "," + str(intersection) + "," + str(union) + "," + str(Total_TRA) + "," + str(TRA_1) + "," + str(TRA_2) + "," + str(TRA_gt2)+ "," + str(Total_TRB) + "," + str(TRB_1) + "," + str(TRB_2) + "," + str(TRB_gt2)+","+str(k_full)+","+str(k_full_match)+","+str(k_N)+","+str(k_N_match))
file.write("\n")





