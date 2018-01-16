#!/bin/bash

source $(dirname $0)/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('barcodes',help='File with database of barcodes. One barcode per line')
parser.add_argument('f2',help='fastq file. R2 file with actual barcodes')
parser.add_argument('imrep_file',help='file obtained with imrep with extended output settings')
parser.add_argument('out',help='Prefix for output')

EOF

dirSource="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "Start analysis"

python ${dirSource}/imrep.10x.step1.py ${BARCODES} ${F2} ${IMREP_FILE}  ${OUT}


#03T3.out_raw_barcode.txt
#03T3.out_raw_barcode2.txt
#03T3.out_raw_barcode2_ed1.txt
#03T3.out_raw_barcode_ed1.txt

python ${dirSource}/imrep.10x.step2.py ${OUT}_raw_barcode.txt ${OUT}.ed0


echo "done!"
