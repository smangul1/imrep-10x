import pysam
import argparse
import csv
import sys
import os
from Bio import SeqIO
from Levenshtein import distance





#Assumes all barcodes in dict_reads_CDR3_TRA mathed R2
#Given dict_reads_CDR3_TRA match to corresponding CDR3s
def check_read_TCR(dict_read_barcode, dict_reads_CDR3_TRA, dict_reads_CDR3_TRB, dict_barcode_TRA_local,dict_barcode_TRB_local,raw_file,flag):
    k=0

    reads_TRA=dict_reads_CDR3_TRA.keys()
    reads_TRB=dict_reads_CDR3_TRB.keys()
    
    print (len(reads_TRA), len(reads_TRB))

    file=open(raw_file,"w")
    
    
    file.write("chain,barcode,read,cdr3,flag")
    file.write("\n")

    for key, value in dict_read_barcode.items():
        read=key
        barcode=value
        if read in reads_TRA:
            if read in dict_reads_CDR3_TRA.keys():
                k+=1
                cdr3=dict_reads_CDR3_TRA[read]
                dict_barcode_TRA_local[barcode].append(cdr3)
                file.write("TRA,"+str(barcode)+","+str(read)+","+str(cdr3)+","+str(flag))
                file.write("\n")
        
     
        elif read in reads_TRB:
            if read in dict_reads_CDR3_TRB.keys():
                k+=1
                cdr3=dict_reads_CDR3_TRB[read]
                dict_barcode_TRB_local[barcode].append(cdr3)
                file.write("TRB,"+str(barcode)+","+str(read)+","+str(cdr3)+","+str(flag))
                file.write("\n")
        

    file.close()
    return k




#search a string from set with ED=1. Report only if there is one such sequence
def search_ed_1 (set_barcodes, cdr3_barcode):
    match=[]
    
    for s in set_barcodes:
        if len(s)==len(cdr3_barcode):
            if distance(str(s),str(cdr3_barcode))==1:
                match.append(s)
        else:
            print ("Warning : length is different",s,cdr3_barcode)



    if len(match)==1:
        return match[0]
    else:
        return ""


#=========================

ap = argparse.ArgumentParser()
ap.add_argument('barcodes', help='File with database of barcodes. One barcode per line')
ap.add_argument('f2', help='fastq file. R2 file with actual barcodes')
ap.add_argument('imrep_file', help='file obtained with imrep with extended output settings')
ap.add_argument('out', help='Prefix for output')


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
    






print ("Number of barcodes", len(bSet), len(bSet2))
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
    
    elif row[2]=="TRB":
         dict_reads_CDR3_TRB[read]=cdr3







print ("Number of TRA reads", len(reads_TRA), len(setCDR3_TRA),sum(dict_cdr3_nReads_TRA.values()))
print ("Number of TRB reads", len(reads_TRB), len(setCDR3_TRB), sum(dict_cdr3_nReads_TRB.values()))









print ("Group CDR3s by barcodes ...")

fileraw1=args.out+"_raw_barcode.txt"
fileraw2=args.out+"_raw_barcode_ed1.txt"
fileraw3=args.out+"_raw_barcode2.txt"
fileraw4=args.out+"_raw_barcode2_ed1.txt"

dict_read_barcode={}

dict_barcode_TRA={}
dict_barcode_TRB={}
dict_barcode_ed1_TRA={}
dict_barcode_ed1_TRB={}
dict_barcode2_TRA={}
dict_barcode2_TRB={}
dict_barcode2_ed1_TRA={}
dict_barcode2_ed1_TRB={}


for b in bSet:
    
    dict_barcode_TRA[b]=[]
    dict_barcode_TRB[b]=[]
    dict_barcode_ed1_TRA[b]=[]
    dict_barcode_ed1_TRB[b]=[]
    dict_barcode2_TRA[b]=[]
    dict_barcode2_TRB[b]=[]
    dict_barcode2_ed1_TRA[b]=[]
    dict_barcode2_ed1_TRB[b]=[]
    




dict_read_barcode={}
dict_read_barcode_ed1={}
dict_read_barcode2={}
dict_read_barcode2_ed1={}


readsSet=set()



with open(args.f2, "rU") as handle:
    for record in SeqIO.parse(handle, "fastq") :
        barcode=str(record.seq)
        read=str(record.id)
        readsSet.add(read)
        
        
        
        
        if barcode[0]=='N': # barcode with N as first letter
            barcodeShort=barcode[1:]
            if barcodeShort in bSet2: # perfect match
                dict_read_barcode2[read]=dictB[barcodeShort]
            
            else:
                correct_barcode=search_ed_1(bSet2,barcodeShort)
                if correct_barcode!="":
                    dict_read_barcode2_ed1[read]=dictB[correct_barcode]
        else: # normal barcode
            if barcode in bSet: # perfect match
                dict_read_barcode[read]=barcode
            else:
                correct_barcode=search_ed_1(bSet,barcode)
                if correct_barcode!="":
                    dict_read_barcode_ed1[read]=correct_barcode






# After we identify to which barcode belong each read, we take those one which have CDR3

k1=0
k2=0
k3=0
k4=0

k1=check_read_TCR(dict_read_barcode, dict_reads_CDR3_TRA, dict_reads_CDR3_TRB, dict_barcode_TRA,dict_barcode_TRB,fileraw1,0)
k2=check_read_TCR(dict_read_barcode_ed1, dict_reads_CDR3_TRA, dict_reads_CDR3_TRB, dict_barcode_ed1_TRA,dict_barcode_ed1_TRB,fileraw2,1)
k3=check_read_TCR(dict_read_barcode2, dict_reads_CDR3_TRA, dict_reads_CDR3_TRB, dict_barcode2_TRA,dict_barcode2_TRB,fileraw3,2)
k4=check_read_TCR(dict_read_barcode2_ed1, dict_reads_CDR3_TRA, dict_reads_CDR3_TRB, dict_barcode2_ed1_TRA,dict_barcode2_ed1_TRB,fileraw4,3)

print (len(readsSet),len(dict_read_barcode)+len(dict_read_barcode_ed1)+ len(dict_read_barcode2)+len(dict_read_barcode2_ed1))
print (len(dict_read_barcode),len(dict_read_barcode_ed1), len(dict_read_barcode2),len(dict_read_barcode2_ed1))

print ("-----")

print (len(dict_read_barcode),k1)
print (len(dict_read_barcode_ed1),k2)
print (len(dict_read_barcode2), k3)
print (len(dict_read_barcode2_ed1), k4)


file=open(args.out+".log","w")

file.write("Total reads, Reads with matchng barcodes full ed=0,Reads with matchng barcodes full ed=1, Reads with matchng barcodes non-full ed=0,Reads with matchng barcodes non-full ed=1,Reads with matchng barcodes full ed=0 TCR,Reads with matchng barcodes full ed=1 TCR, Reads with matchng barcodes non-full ed=0 TCR,Reads with matchng barcodes non-full ed=1 TCR")
file.write("\n")


file.write(str(len(readsSet))+","+str(len(dict_read_barcode))+","+str(len(dict_read_barcode_ed1))+","+str(len(dict_read_barcode2))+","+str(len(dict_read_barcode2_ed1))+","+str(k1)+","+str(k2)+","+str(k3)+","+str(k4)   )



file.close()



