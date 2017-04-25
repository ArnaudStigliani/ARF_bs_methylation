#!/bin/bash

intersectBed -a mCG.bedGraph -b G2_ARF2.bed > G2_CG_met_ARF2.bed &
intersectBed -a mCHG.bedGraph -b G2_ARF2.bed > G2_CHG_met_ARF2.bed &
intersectBed -a mCHH.bedGraph -b G2_ARF2.bed > G2_CHH_met_ARF2.bed &
wait
intersectBed -a mCG.bedGraph -b C4_ARF2.bed > C4_CG_met_ARF2.bed &
intersectBed -a mCHG.bedGraph -b C4_ARF2.bed > C4_CHG_met_ARF2.bed &
intersectBed -a mCHH.bedGraph -b C4_ARF2.bed > C4_CHH_met_ARF2.bed &
wait
intersectBed -a mCG.bedGraph -b G5_ARF2.bed > G5_CG_met_ARF2.bed &
intersectBed -a mCHG.bedGraph -b G5_ARF2.bed > G5_CHG_met_ARF2.bed &
intersectBed -a mCHH.bedGraph -b G5_ARF2.bed > G5_CHH_met_ARF2.bed &
wait
intersectBed -a mCG.bedGraph -b G6_ARF2.bed > G6_CG_met_ARF2.bed &
intersectBed -a mCHG.bedGraph -b G6_ARF2.bed > G6_CHG_met_ARF2.bed &
intersectBed -a mCHH.bedGraph -b G6_ARF2.bed > G6_CHH_met_ARF2.bed &

exit 0
