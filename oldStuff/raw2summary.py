import csv
import argparse
import sys
import os

#=========================

ap = argparse.ArgumentParser()
ap.add_argument('raw', help='raw file produced by barcode.py')
ap.add_argument('out', help='raw file produced by barcode.py')

args = ap.parse_args()


#chain,barcode,read,cdr3,flag
#TRB,TGACTTTGGACGGA,M02564:19:000000000-G15E1:1:1103:17105:15203,CASSPYGNTIYF,0


base=os.path.basename(args.raw)
sample=os.path.splitext(base)[0]



file=open(args.raw)
csv_f=csv.reader(file)
next(csv_f, None)



barcode_dict={}
barcode_set=set()


for row in csv_f:
    chain=row[0]
    barcode=row[1]
    cdr3=row[3]
    
    
    
    
    if barcode in barcode_set:
        if chain=="TRA":
            barcode_dict[barcode][0].append(cdr3)
        else:
            barcode_dict[barcode][1].append(cdr3)
    else:
        barcode_set.add(barcode)
        barcode_dict[barcode]=([],[])
        if chain=="TRA":
            barcode_dict[barcode][0].append(cdr3)
        else:
            barcode_dict[barcode][1].append(cdr3)





#print set(barcode_TRB_dict['TGACTTTGGACGGA'])
#print set(barcode_TRB_dict['TTGGTACTTTACTC'])

#write output --------------------------------------------

file=open(args.out,"w")


for key, value in barcode_dict.iteritems():
    
    if len(set(value[0]))>1 or len(set(value[1]))>1:
        print len(set(value[0])),set(value[0]),set(value[1])
    
    
        print sample,key, "TCRA",';'.join(list(set(value[0])))
        print sample,key, "TCRB",';'.join(list(set(value[1])))
    
        for s in set(value[0]):
            print "---",value[0].count(s), s
    
    else:
        TCRA_cdr3=list(set(value[0]))
        TCRB_cdr3=list(set(value[1]))
        
        if len(TCRA_cdr3)>0 and len(TCRB_cdr3)>0:
            print sample,key,TCRA_cdr3[0],TCRB_cdr3[0], len(value[0]), len(value[1])
        elif len(TCRA_cdr3)>0:
            print sample,key,TCRA_cdr3[0],"N/A", len(value[0]), "N/A"
        elif len(TCRB_cdr3)>0:
            print sample,key,"N/A",TCRB_cdr3[0], "N/A", len(value[1])





#file.write(key+",0,1,1"+)


