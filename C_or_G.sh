#!/bin/bash

awk 'BEGIN{OFS="\t"}{print $5,$6,$7,$8}' G2_C_all_met_ARF5_pos.bed | sed 's/-1$/1\t1\t-/' | sed 's/1$/1\t1\t+/'| bedtools getfasta -fi ~/Data/tair10/tair10.fas -fo lol.fas -bed -  -s

exit 0
