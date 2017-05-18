#!/bin/bash

bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed neg_Ind1_neg.bed -fo neg_Ind1_neg.fas
bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed neg_Ind2_neg.bed -fo neg_Ind2_neg.fas
bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed neg_Ind3_neg.bed -fo neg_Ind3_neg.fas
bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed neg_Ind4_neg.bed -fo neg_Ind4_neg.fas
bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed promoters_Aux_30_Ind.bed -fo promoters_Aux_30_Ind.fas
bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed neg_Rep1_neg.bed -fo neg_Rep1_neg.fas
bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed neg_Rep2_neg.bed -fo neg_Rep2_neg.fas
bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed neg_Rep3_neg.bed -fo neg_Rep3_neg.fas
bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed neg_Rep4_neg.bed -fo neg_Rep4_neg.fas
bedtools getfasta -fi ~/Data/tair10/tair10.fas -bed promoters_Aux_30_Rep.bed -fo promoters_Aux_30_Rep.fas

exit 0
