#!/bin/bash

length=$1

./get_promoter.sh $length
python ~/negative_arnaud/fonctions.py --pos promoters_Aux_30_Ind.bed --neg neg_Ind -n 4
python ~/negative_arnaud/fonctions.py --pos promoters_Aux_30_Rep.bed --neg neg_Rep -n 4
./getfasta.sh
Rscript scores.r Ind
Rscript scores.r Rep

exit 0
