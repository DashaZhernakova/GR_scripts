#!/bin/bash

f=$1
r=0.2
if [ "$#" -gt 1 ]; then
 r=$2
fi

echo "Pruning using 1000 50 ${r}"
plink \
--bfile $f \
--indep-pairwise 1000 50 ${r} \
--out ${f}.pruned_r2_${r}

plink \
--bfile $f \
--extract ${f}.pruned_r2_${r}.prune.in \
--make-bed \
--out ${f}.pruned_r2_${r}
