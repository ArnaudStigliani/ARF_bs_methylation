#!/bin/bash


intersectBed -a mCG.bedGraph -b T1_ARF2_pos.bed > T1_CG_met_ARF2_pos.bed &
intersectBed -a mCHG.bedGraph -b T1_ARF2_pos.bed > T1_CHG_met_ARF2_pos.bed &
intersectBed -a mCHH.bedGraph -b T1_ARF2_pos.bed > T1_CHH_met_ARF2_pos.bed &
#
intersectBed -a mCG.bedGraph -b G2_ARF2_pos.bed > G2_CG_met_ARF2_pos.bed &
intersectBed -a mCHG.bedGraph -b G2_ARF2_pos.bed > G2_CHG_met_ARF2_pos.bed &
intersectBed -a mCHH.bedGraph -b G2_ARF2_pos.bed > G2_CHH_met_ARF2_pos.bed &
#
intersectBed -a mCG.bedGraph -b T3_ARF2_pos.bed > T3_CG_met_ARF2_pos.bed &
intersectBed -a mCHG.bedGraph -b T3_ARF2_pos.bed > T3_CHG_met_ARF2_pos.bed &
intersectBed -a mCHH.bedGraph -b T3_ARF2_pos.bed > T3_CHH_met_ARF2_pos.bed &
#
intersectBed -a mCG.bedGraph -b C4_ARF2_pos.bed > C4_CG_met_ARF2_pos.bed &
intersectBed -a mCHG.bedGraph -b C4_ARF2_pos.bed > C4_CHG_met_ARF2_pos.bed &
intersectBed -a mCHH.bedGraph -b C4_ARF2_pos.bed > C4_CHH_met_ARF2_pos.bed &
#
intersectBed -a mCG.bedGraph -b G5_ARF2_pos.bed > G5_CG_met_ARF2_pos.bed &
intersectBed -a mCHG.bedGraph -b G5_ARF2_pos.bed > G5_CHG_met_ARF2_pos.bed &
intersectBed -a mCHH.bedGraph -b G5_ARF2_pos.bed > G5_CHH_met_ARF2_pos.bed &
#
intersectBed -a mCG.bedGraph -b G6_ARF2_pos.bed > G6_CG_met_ARF2_pos.bed &
intersectBed -a mCHG.bedGraph -b G6_ARF2_pos.bed > G6_CHG_met_ARF2_pos.bed &
intersectBed -a mCHH.bedGraph -b G6_ARF2_pos.bed > G6_CHH_met_ARF2_pos.bed &
wait
intersectBed -a mCG.bedGraph -b T1_ARF2_neg.bed > T1_CG_met_ARF2_neg.bed &
intersectBed -a mCHG.bedGraph -b T1_ARF2_neg.bed > T1_CHG_met_ARF2_neg.bed &
intersectBed -a mCHH.bedGraph -b T1_ARF2_neg.bed > T1_CHH_met_ARF2_neg.bed &
#
intersectBed -a mCG.bedGraph -b G2_ARF2_neg.bed > G2_CG_met_ARF2_neg.bed &
intersectBed -a mCHG.bedGraph -b G2_ARF2_neg.bed > G2_CHG_met_ARF2_neg.bed &
intersectBed -a mCHH.bedGraph -b G2_ARF2_neg.bed > G2_CHH_met_ARF2_neg.bed &
#
intersectBed -a mCG.bedGraph -b T3_ARF2_neg.bed > T3_CG_met_ARF2_neg.bed &
intersectBed -a mCHG.bedGraph -b T3_ARF2_neg.bed > T3_CHG_met_ARF2_neg.bed &
intersectBed -a mCHH.bedGraph -b T3_ARF2_neg.bed > T3_CHH_met_ARF2_neg.bed &
#
intersectBed -a mCG.bedGraph -b C4_ARF2_neg.bed > C4_CG_met_ARF2_neg.bed &
intersectBed -a mCHG.bedGraph -b C4_ARF2_neg.bed > C4_CHG_met_ARF2_neg.bed &
intersectBed -a mCHH.bedGraph -b C4_ARF2_neg.bed > C4_CHH_met_ARF2_neg.bed &
#
intersectBed -a mCG.bedGraph -b G5_ARF2_neg.bed > G5_CG_met_ARF2_neg.bed &
intersectBed -a mCHG.bedGraph -b G5_ARF2_neg.bed > G5_CHG_met_ARF2_neg.bed &
intersectBed -a mCHH.bedGraph -b G5_ARF2_neg.bed > G5_CHH_met_ARF2_neg.bed &
#
intersectBed -a mCG.bedGraph -b G6_ARF2_neg.bed > G6_CG_met_ARF2_neg.bed &
intersectBed -a mCHG.bedGraph -b G6_ARF2_neg.bed > G6_CHG_met_ARF2_neg.bed &
intersectBed -a mCHH.bedGraph -b G6_ARF2_neg.bed > G6_CHH_met_ARF2_neg.bed &
wait
intersectBed -a mCG.bedGraph -b T1_ARF5_pos.bed > T1_CG_met_ARF5_pos.bed &
intersectBed -a mCHG.bedGraph -b T1_ARF5_pos.bed > T1_CHG_met_ARF5_pos.bed &
intersectBed -a mCHH.bedGraph -b T1_ARF5_pos.bed > T1_CHH_met_ARF5_pos.bed &
#
intersectBed -a mCG.bedGraph -b G2_ARF5_pos.bed > G2_CG_met_ARF5_pos.bed &
intersectBed -a mCHG.bedGraph -b G2_ARF5_pos.bed > G2_CHG_met_ARF5_pos.bed &
intersectBed -a mCHH.bedGraph -b G2_ARF5_pos.bed > G2_CHH_met_ARF5_pos.bed &
#
intersectBed -a mCG.bedGraph -b T3_ARF5_pos.bed > T3_CG_met_ARF5_pos.bed &
intersectBed -a mCHG.bedGraph -b T3_ARF5_pos.bed > T3_CHG_met_ARF5_pos.bed &
intersectBed -a mCHH.bedGraph -b T3_ARF5_pos.bed > T3_CHH_met_ARF5_pos.bed &
#
intersectBed -a mCG.bedGraph -b C4_ARF5_pos.bed > C4_CG_met_ARF5_pos.bed &
intersectBed -a mCHG.bedGraph -b C4_ARF5_pos.bed > C4_CHG_met_ARF5_pos.bed &
intersectBed -a mCHH.bedGraph -b C4_ARF5_pos.bed > C4_CHH_met_ARF5_pos.bed &
#
intersectBed -a mCG.bedGraph -b G5_ARF5_pos.bed > G5_CG_met_ARF5_pos.bed &
intersectBed -a mCHG.bedGraph -b G5_ARF5_pos.bed > G5_CHG_met_ARF5_pos.bed &
intersectBed -a mCHH.bedGraph -b G5_ARF5_pos.bed > G5_CHH_met_ARF5_pos.bed &
#
intersectBed -a mCG.bedGraph -b G6_ARF5_pos.bed > G6_CG_met_ARF5_pos.bed &
intersectBed -a mCHG.bedGraph -b G6_ARF5_pos.bed > G6_CHG_met_ARF5_pos.bed &
intersectBed -a mCHH.bedGraph -b G6_ARF5_pos.bed > G6_CHH_met_ARF5_pos.bed &
wait
intersectBed -a mCG.bedGraph -b T1_ARF5_neg.bed > T1_CG_met_ARF5_neg.bed &
intersectBed -a mCHG.bedGraph -b T1_ARF5_neg.bed > T1_CHG_met_ARF5_neg.bed &
intersectBed -a mCHH.bedGraph -b T1_ARF5_neg.bed > T1_CHH_met_ARF5_neg.bed &
#
intersectBed -a mCG.bedGraph -b G2_ARF5_neg.bed > G2_CG_met_ARF5_neg.bed &
intersectBed -a mCHG.bedGraph -b G2_ARF5_neg.bed > G2_CHG_met_ARF5_neg.bed &
intersectBed -a mCHH.bedGraph -b G2_ARF5_neg.bed > G2_CHH_met_ARF5_neg.bed &
#
intersectBed -a mCG.bedGraph -b T3_ARF5_neg.bed > T3_CG_met_ARF5_neg.bed &
intersectBed -a mCHG.bedGraph -b T3_ARF5_neg.bed > T3_CHG_met_ARF5_neg.bed &
intersectBed -a mCHH.bedGraph -b T3_ARF5_neg.bed > T3_CHH_met_ARF5_neg.bed &
#
intersectBed -a mCG.bedGraph -b C4_ARF5_neg.bed > C4_CG_met_ARF5_neg.bed &
intersectBed -a mCHG.bedGraph -b C4_ARF5_neg.bed > C4_CHG_met_ARF5_neg.bed &
intersectBed -a mCHH.bedGraph -b C4_ARF5_neg.bed > C4_CHH_met_ARF5_neg.bed &
#
intersectBed -a mCG.bedGraph -b G5_ARF5_neg.bed > G5_CG_met_ARF5_neg.bed &
intersectBed -a mCHG.bedGraph -b G5_ARF5_neg.bed > G5_CHG_met_ARF5_neg.bed &
intersectBed -a mCHH.bedGraph -b G5_ARF5_neg.bed > G5_CHH_met_ARF5_neg.bed &
#
intersectBed -a mCG.bedGraph -b G6_ARF5_neg.bed > G6_CG_met_ARF5_neg.bed &
intersectBed -a mCHG.bedGraph -b G6_ARF5_neg.bed > G6_CHG_met_ARF5_neg.bed &
intersectBed -a mCHH.bedGraph -b G6_ARF5_neg.bed > G6_CHH_met_ARF5_neg.bed &
wait

exit 0
