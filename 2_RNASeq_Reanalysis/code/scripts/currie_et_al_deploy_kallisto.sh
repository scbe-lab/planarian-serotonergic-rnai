#!/bin/bash
KALLISTO_INDEX="/mnt/sda/alberto/colabos/elena_sero/rnaseq_reanalysis/data/genome/Smed.fna.kallisto"
KALLISTO_OUTPUT_FOLDER="/mnt/sda/alberto/colabos/elena_sero/rnaseq_reanalysis/outputs/RNA/output_currie/kallisto_output/raw_output/"
SAMPLE_FOLDER="/mnt/sda/alberto/colabos/elena_sero/rnaseq_reanalysis/data/currie_et_al_2013/fq"

for f in ${SAMPLE_FOLDER}/* ; do
    x=${f##*/}
    echo Starting with sample ${x} ...
    mkdir -p ${KALLISTO_OUTPUT_FOLDER}/kallisto_out_${x}
    kallisto quant -t 10 -i \
    $KALLISTO_INDEX \
    --single $f/*.f*.gz \
    -l 200 -s 30 \
    -o ${KALLISTO_OUTPUT_FOLDER}/kallisto_out_${x}
    echo "done with sample ${x} ..."
done
echo "Done."
