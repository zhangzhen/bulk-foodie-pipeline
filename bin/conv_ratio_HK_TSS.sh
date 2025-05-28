#!/bin/bash

# human
TSS_HK_1k=/dshare/home/xiec/Research/prj/brain/rsl/5409/TSS_HK/ENCFF493CCB_HK_1k.bed
#TSS_HK_1k=/lustre/grp/xxllab/xiec/Research/prj/brain/rsl/5409/TSS_HK/ENCFF493CCB_HK_1k.bed
# mouse
#TSS_HK_1k=/dshare/home/xiec/Research/prj/brain/rsl/5409/TSS_HK/ENCFF498BEJ_HK_1k.bed
#TSS_HK_1k=/lustre/grp/xxllab/xiec/Research/prj/brain/rsl/5409/TSS_HK/ENCFF498BEJ_HK_1k.bed

tmp=tmp_TSS
mkdir -p ${tmp}/

IN=../02-merge-outputs/ratio
OUT=../02-merge-outputs/ratio_TSS

mkdir -p ${OUT}/

ls $IN | while read f; do
    mkdir -p ${OUT}/${f}/
    echo "#!/bin/bash" > ${tmp}/${f}.sh
    echo "bedtools intersect -a ${IN}/${f}/${f}_sites.txt -b $TSS_HK_1k -wa -wb | awk 'BEGIN {OFS=\"\\t\"} {if (\$13 == \"+\") rel_pos = \$3 - \$10 + 500; else rel_pos = \$10 - \$3 - 500; print \$4, \$5, rel_pos}' > ${OUT}/${f}/${f}_TSS_HK.tsv" >> ${tmp}/${f}.sh
    sbatch -t 10:00:00 -p compute_new -c 2 --mem=20g -e ${tmp}/${f}.e -o ${tmp}/${f}.o ${tmp}/${f}.sh
done

