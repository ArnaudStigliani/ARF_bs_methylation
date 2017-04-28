#!/bin/bash


intersectBed -wb -a mC_all.bedGraph -b T1_ARF2_pos.bed > T1_C_all_met_ARF2_pos.bed &
intersectBed -wb -a mC_all.bedGraph -b G2_ARF2_pos.bed > G2_C_all_met_ARF2_pos.bed &
intersectBed -wb -a mC_all.bedGraph -b T3_ARF2_pos.bed > T3_C_all_met_ARF2_pos.bed &
intersectBed -wb -a mC_all.bedGraph -b C4_ARF2_pos.bed > C4_C_all_met_ARF2_pos.bed &
intersectBed -wb -a mC_all.bedGraph -b G5_ARF2_pos.bed > G5_C_all_met_ARF2_pos.bed &
intersectBed -wb -a mC_all.bedGraph -b G6_ARF2_pos.bed > G6_C_all_met_ARF2_pos.bed &
intersectBed -wb -a mC_all.bedGraph -b T1_ARF2_neg.bed > T1_C_all_met_ARF2_neg.bed &
intersectBed -wb -a mC_all.bedGraph -b G2_ARF2_neg.bed > G2_C_all_met_ARF2_neg.bed &
intersectBed -wb -a mC_all.bedGraph -b T3_ARF2_neg.bed > T3_C_all_met_ARF2_neg.bed &
intersectBed -wb -a mC_all.bedGraph -b C4_ARF2_neg.bed > C4_C_all_met_ARF2_neg.bed &
intersectBed -wb -a mC_all.bedGraph -b G5_ARF2_neg.bed > G5_C_all_met_ARF2_neg.bed &
intersectBed -wb -a mC_all.bedGraph -b G6_ARF2_neg.bed > G6_C_all_met_ARF2_neg.bed &
intersectBed -wb -a mC_all.bedGraph -b T1_ARF5_pos.bed > T1_C_all_met_ARF5_pos.bed &
intersectBed -wb -a mC_all.bedGraph -b G2_ARF5_pos.bed > G2_C_all_met_ARF5_pos.bed &
intersectBed -wb -a mC_all.bedGraph -b T3_ARF5_pos.bed > T3_C_all_met_ARF5_pos.bed &
intersectBed -wb -a mC_all.bedGraph -b C4_ARF5_pos.bed > C4_C_all_met_ARF5_pos.bed &
intersectBed -wb -a mC_all.bedGraph -b G5_ARF5_pos.bed > G5_C_all_met_ARF5_pos.bed &
intersectBed -wb -a mC_all.bedGraph -b G6_ARF5_pos.bed > G6_C_all_met_ARF5_pos.bed &
intersectBed -wb -a mC_all.bedGraph -b T1_ARF5_neg.bed > T1_C_all_met_ARF5_neg.bed &
intersectBed -wb -a mC_all.bedGraph -b G2_ARF5_neg.bed > G2_C_all_met_ARF5_neg.bed &
intersectBed -wb -a mC_all.bedGraph -b T3_ARF5_neg.bed > T3_C_all_met_ARF5_neg.bed &
intersectBed -wb -a mC_all.bedGraph -b C4_ARF5_neg.bed > C4_C_all_met_ARF5_neg.bed &
intersectBed -wb -a mC_all.bedGraph -b G5_ARF5_neg.bed > G5_C_all_met_ARF5_neg.bed &
intersectBed -wb -a mC_all.bedGraph -b G6_ARF5_neg.bed > G6_C_all_met_ARF5_neg.bed &
wait
sort -k1,1 -k2,2n sites_ARF5_neg.bed | mergeBed -i - | intersectBed -v -a ARF5_neg.bed -b - | intersectBed -a mC_all.bedGraph -b - > reg_wo_sites_C_all_met_ARF5_neg.bed &
sort -k1,1 -k2,2n sites_ARF5_pos.bed | mergeBed -i - | intersectBed -v -a ARF5_pos.bed -b - | intersectBed -a mC_all.bedGraph -b - > reg_wo_sites_C_all_met_ARF5_pos.bed &
sort -k1,1 -k2,2n sites_ARF2_neg.bed | mergeBed -i - | intersectBed -v -a ARF2_neg.bed -b - | intersectBed -a mC_all.bedGraph -b - > reg_wo_sites_C_all_met_ARF2_neg.bed &
sort -k1,1 -k2,2n sites_ARF2_pos.bed | mergeBed -i - | intersectBed -v -a ARF2_pos.bed -b - | intersectBed -a mC_all.bedGraph -b - > reg_wo_sites_C_all_met_ARF2_pos.bed &
exit 0
