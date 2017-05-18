#!/bin/bash

length=$1

awk -v var=$length 'BEGIN{FS=";";OFS="\t"} ($5~"+") {print ($2,$3-var,$3)} ($5~"-") {print ($2,$4,$4+var)}' Aux_30_Ind_Protein_Models.csv | sed -e 's/Chr/chr/' -e '/chr[1-5]/!d'  > promoters_Aux_30_Ind.bed

awk -v var=$length 'BEGIN{FS=";";OFS="\t"} ($5~"+") {print ($2,$3-var,$3)} ($5~"-") {print ($2,$4,$4+var)}' Aux_30_Rep_Proteins_Models.csv | sed -e 's/Chr/chr/' -e '/chr[1-5]/!d' > promoters_Aux_30_Rep.bed

exit 0
