import pysam
import argparse
import csv
import sys
import os
from Bio import SeqIO
from Levenshtein import distance





#Assumes all barcodes in dict_reads_CDR3_TRA mathed R2
#Given dict_reads_CDR3_TRA match to corresponding CDR3s


def check_read_TCR(barcodesSet,dict_read_barcode, reads_TRA, reads_TRB,flag, dict_reads_CDR3_TRA, dict_reads_CDR3_TRB, dict_barcode_TRA,dict_barcode_TRB,raw_file):

    file=open(raw_file,"w")

    for key, value in dict_read_barcode.iteritems():
        if flag==0:
            barcode=value[1:]
        else:
            barcode=value
        
        read=key

    if read in reads_TRA:
        file.write("TRA,"+str(barcode)+","+dictB[barcode]+","+str(read)+","+str(cdr3)+","+ str(flag))
        file.write("\n")
        cdr3=dict_reads_CDR3_TRA[read]
        dict_barcode_TRA[barcodeShort].append(cdr3)
 
    elif read in reads_TRB:
        fRaw2.write("TRB,"+str(barcodeShort)+","+dictB[barcodeShort]+","+str(read)+","+ str(flag))
        fRaw2.write("\n")
        cdr3=dict_reads_CDR3_TRB[read]
        dict_barcode_TRB[barcodeShort].append(cdr3)






#search a string from set with ED=1. Report only if there is one such sequence
def search_ed_1 (set_barcodes, cdr3_barcode):
    match=[]
    
    for s in set_barcodes:
        if len(s)==len(cdr3_barcode):
            if distance(str(s),str(cdr3_barcode))==1:
                match.append(cdr3_barcode)
        else:
            print "Warning : length is different",s,cdr3_barcode



    if len(match)==1:
        return match[0]
    else:
        return False




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
dictB={}

for row in csv_f:
    b=row[0].split("-")[0]
    bSet.add(b)
    bSet2.add(b[1:])
    dictB[b[1:]]=b
    






print "Number of barcodes", len(bSet), len(bSet2)
k_full=0
k_full_match=0
k_full_match_ed1=0
k_full_match_TCR=0
k_full_match_TCR_ed1=0
k_N=0
k_N_match=0
k_N_match_ed1=0
k_N_match_TCR=0
k_N_match_TCR_ed1=0


barcodesSet=set()

#CDR3
#M02564:19:000000000-G15E1:1:2102:15042:11150    CILGSNNNDMRF    TRA     TRAV26  NA      TRAV26-2*01:20:2        TRAJ43*01:17:0  1       1       1


dict_reads_CDR3_TRA={}
dict_reads_CDR3_TRB={}


setCDR3_TRA=set()
setCDR3_TRB=set()

dict_cdr3_nReads_TRA={}
dict_cdr3_nReads_TRB={}


reads_TRA=set()
reads_TRB=set()




#imrep file -----------------------------------


n_CDR3_imrep=0

f=open(args.imrep_file,"r")
csv_f=csv.reader(f,delimiter="\t")
next(csv_f, None)  # skip the headers
    

for row in csv_f:
    n_CDR3_imrep+=1
    read=row[0]
    cdr3=row[1]
    if row[2]=="TRA":
        dict_reads_CDR3_TRA[read]=cdr3
        reads_TRA.add(read)
        
        if cdr3 not in setCDR3_TRA:
            dict_cdr3_nReads_TRA[cdr3]=1
        else:
            dict_cdr3_nReads_TRA[cdr3]+=1
        
        setCDR3_TRA.add(cdr3)
    elif row[2]=="TRB":
         dict_reads_CDR3_TRB[read]=cdr3
         reads_TRB.add(read)
         
         if cdr3 not in setCDR3_TRB:
             dict_cdr3_nReads_TRB[cdr3]=1
         else:
             dict_cdr3_nReads_TRB[cdr3]+=1
         setCDR3_TRB.add(cdr3)

    






print "Number of TRA reads", len(reads_TRA), len(setCDR3_TRA),sum(dict_cdr3_nReads_TRA.values())
print "Number of TRB reads", len(reads_TRB), len(setCDR3_TRB), sum(dict_cdr3_nReads_TRB.values())



fileraw2=args.out+"_raw2.txt"





print "Group CDR3s by barcodes ..."

dict_read_barcode={}

dict_barcode_TRA={}
dict_barcode_TRB={}
dict_barcode2_TRA={}
dict_barcode2_TRB={}


dict_read_barcode={}
dict_read_barcode2={}
with open(args.f2, "rU") as handle:
    for record in SeqIO.parse(handle, "fastq") :
        barcode=record.seq
        read=record.id
        
        if barcode[0]=='N':
            print barcode
        


sys.exit(1)


#bSet - full barcodes
for b in bSet:
    dict_barcode_TRA[b]=[]
    dict_barcode_TRB[b]=[]

#bSet - short barcodes i.e. excluding first character
for b in bSet2:
    dict_barcode2_TRA[b]=[]
    dict_barcode2_TRB[b]=[]




#check_read_TCR(barcodesSet,dict_read_barcode, reads_TRA, reads_TRB,flag, dict_reads_CDR3_TRA, dict_reads_CDR3_TRB, dict_barcode_TRA,dict_barcode_TRB raw_file):




with open(args.f2, "rU") as handle:
    for record in SeqIO.parse(handle, "fastq") :
        
        barcodeShort=record.seq[1:]
        barcode=record.seq


        read=record.id
        

        if barcodeShort in bSet2:
                    barcodesSet.add(barcodeShort)
                    k_N_match+=1
                    if read in reads_TRA:
                        k_N_match_TCR+=1
                        fRaw2.write("TRA,"+str(barcodeShort)+","+dictB[barcodeShort]+","+str(read)+","+str(cdr3)+",2")
                        fRaw2.write("\n")
                        cdr3=dict_reads_CDR3_TRA[read]
                        dict_barcode_TRA[barcodeShort].append(cdr3)
                        set_barcode_TRA.add(barcodeShort)
                            #print barcodeShort,cdr3,read
                            #sys.exit(1)
                    elif read in reads_TRB:
                        k_N_match_TCR+=1
                        fRaw2.write("TRB,"+str(barcodeShort)+","+dictB[barcodeShort]+","+str(read)+","+str(cdr3)+",2")
                        fRaw2.write("\n")


                        cdr3=dict_reads_CDR3_TRB[read]
                        dict_barcode_TRB[barcodeShort].append(cdr3)
                        set_barcode_TRB.add(barcodeShort)
                # no perfect match
                elif search_ed_1(bSet2,barcodeShort):
                        barcodesSet.add(barcodeShort)
                        k_N_match_ed1+=1
                        if read in reads_TRA:
                            k_N_match_TCR_ed1+=1
                            fRaw2.write("TRB,"+str(barcodeShort)+","+dictB[barcodeShort]+","+str(read)+","+str(cdr3)+",3")
                            fRaw2.write("\n")
                            cdr3=dict_reads_CDR3_TRA[read]
                            dict_barcode_TRA[barcodeShort].append(cdr3)
                            set_barcode_TRA.add(barcodeShort)
                        elif read in reads_TRB:
                            k_N_match_TCR_ed1+=1
                            fRaw2.write("TRB,"+str(barcodeShort)+","+dictB[barcodeShort]+","+str(read)+","+str(cdr3)+",3")
                            fRaw2.write("\n")
                            cdr3=dict_reads_CDR3_TRB[read]
                            dict_barcode_TRB[barcodeShort].append(cdr3)
                            set_barcode_TRB.add(barcodeShort)
                            
                        
                        
                        
                    





        
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
                    k_full_match_TCR+=1
                    fRaw2.write("TRA,"+str(barcode)+","+str(barcode)+","+str(read)+","+str(cdr3)+",1")
                    fRaw2.write("\n")

                    cdr3=dict_reads_CDR3_TRA[read]
                    dict_barcode_TRA[barcodeShort].append(cdr3)
                    set_barcode_TRA.add(barcodeShort)

                elif read in reads_TRB:
                    k_full_match_TCR+=1
                    fRaw2.write("TRB,"+str(barcode)+","+str(barcode)+","+str(read)+","+str(cdr3)+",1")
                    fRaw2.write("\n")

                    cdr3=dict_reads_CDR3_TRB[read]
                    dict_barcode_TRB[barcodeShort].append(cdr3)
                    set_barcode_TRB.add(barcodeShort)
            elif search_ed_1(bSet,barcode):
                barcodesSet.add(barcode[1:])
                k_full_match_ed1+=1
                if read in reads_TRA:
                    k_full_match_TCR_ed1++1
                    fRaw2.write("TRA,"+str(barcode)+","+str(barcode)+","+str(read)+","+str(cdr3)+",1")
                    fRaw2.write("\n")
                    cdr3=dict_reads_CDR3_TRA[read]
                    dict_barcode_TRA[barcodeShort].append(cdr3)
                    set_barcode_TRA.add(barcodeShort)
            
                elif read in reads_TRB:
                    k_full_match_TCR_ed1+=1
                    fRaw2.write("TRB,"+str(barcode)+","+str(barcode)+","+str(read)+","+str(cdr3)+",1")
                    fRaw2.write("\n")
                    cdr3=dict_reads_CDR3_TRB[read]
                    dict_barcode_TRB[barcodeShort].append(cdr3)
                    set_barcode_TRB.add(barcodeShort)








fRaw2.close()


print "Number of reads with full barcode",k_full
print "Number of reads with full barcode mathing",k_full_match, k_full_match_TCR
print "Number of reads with not-full barcode",k_N
print "Number of reads with not-full barcode mathing",k_N_match, k_N_match_TCR




nCDR3_TRA=[]
nCDR3_TRB=[]

file=open(args.out+"_raw.txt","w")

for key, value in dict_barcode_TRA.iteritems():
    if len(value)>0:
        for s in set(value):
            file.write(base2 +","+dictB[key]+","+"TCRA,"+key +","+s+","+str(len(set(value)))+ ","+str(dict_cdr3_nReads_TRA[s]) )
            file.write("\n")
        nCDR3_TRA.append(len(set(value)))



for key, value in dict_barcode_TRB.iteritems():
    if len(value)>0:
        for s in set(value):
            file.write(base2 +","+dictB[key]+","+"TCRB,"+key +","+s+","+str(len(set(value))) + ","+str(dict_cdr3_nReads_TRB[s]))

            file.write("\n")
        nCDR3_TRB.append(len(set(value)))

file.close()



print "TRA",sum(dict_cdr3_nReads_TRA.values())
print "TRB",sum(dict_cdr3_nReads_TRB.values())




#TRA
file=open(args.out+"_TCR_summary.txt","w")
file.write("Sample, total_barcodes_library, total_barcodes_library2, total_barcodes_R2, intersection_barcodes, union_barcodes,  Total_TRA, TRA_1, TRA_2, TRA>=2,Total_TRB, TRB_1, TRB_2, TRB>=2, Number of reads with full barcode, Number of reads with full barcode mathing, Number of reads with not-full barcode, Number of reads with not-full barcode TCR, Number of reads with not-full barcode mathing, Number of reads with not-full barcode mathing TCR")
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



file.write(base2 +"," + str(len(bSet)) +","+ str(len(bSet2))+ ","+str(total_barcodes) +  "," + str(intersection) + "," + str(union) + "," + str(Total_TRA) + "," + str(TRA_1) + "," + str(TRA_2) + "," + str(TRA_gt2)+ "," + str(Total_TRB) + "," + str(TRB_1) + "," + str(TRB_2) + "," + str(TRB_gt2)+","+str(k_full)+","+str(k_full_match)+","+str(k_full_match_TCR) + ","+str(k_N)+","+str(k_N_match)+","+str(k_N_match_TCR))
file.write("\n")


print n_CDR3_imrep,k_full_match_TCR+k_N_match_TCR


