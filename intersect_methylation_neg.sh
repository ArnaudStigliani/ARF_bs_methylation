#!/bin/bash

intersectBed -a mCG.bedGraph -b ARF2_neg.bed > CG_met_neg_ARF2.bed &
intersectBed -a mCHG.bedGraph -b ARF2_neg.bed > CHG_met_neg_ARF2.bed &
intersectBed -a mCHH.bedGraph -b ARF2_neg.bed > CHH_met_neg_ARF2.bed &

exit 0
