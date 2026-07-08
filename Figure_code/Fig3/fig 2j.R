

#!/usr/bin/env bash
# Usage: bash run_ggsashimi_age.sh <locus> <gene> <celltype> <sample_reads>

LOCUS="$1"
GENE="$2"
CELLTYPE="$3"
SAMPLE_READS="$4"

BASE_DIR="bam"
BAM_DIR="sodu"
IND_AGE="${BASE_DIR}/summary_age.txt_${CELLTYPE}"
OUT_DIR="${BASE_DIR}/${GENE}"
INDIV_DIR="${OUT_DIR}/all_individuals"
OUTPUT_TSV="input_bams.${CELLTYPE}.${GENE}.tsv"
TMP_TSV="${OUTPUT_TSV}.tmp"

mkdir -p "${INDIV_DIR}"
> "${TMP_TSV}"

echo "========================================"
echo "  Locus     : ${LOCUS}"
echo "  Gene      : ${GENE}"
echo "  Cell type : ${CELLTYPE}"
echo "  Reads     : ${SAMPLE_READS}"
echo "========================================"

for AGE in $(cut -d' ' -f2 "${IND_AGE}" | sort -u); do

    echo ""
    echo "--- Age group: ${AGE} ---"

    # Step 1: extract reads per individual
    for IND in $(grep "${AGE}" "${IND_AGE}" | cut -d' ' -f1); do
        INPUT_BAM="${BAM_DIR}/${IND}@${CELLTYPE}.${IND}@${CELLTYPE}.bam"
        OUTPUT_BAM="${INDIV_DIR}/${IND}.${CELLTYPE}.${AGE}.${GENE}.bam"
        echo "  -> ${IND}"
        samtools view "${INPUT_BAM}" "${LOCUS}" -o "${OUTPUT_BAM}"
    done

    # Step 2: merge all individual BAMs
    MERGED_BAM="${OUT_DIR}/mergeAll.${CELLTYPE}.${AGE}.${GENE}.bam"
    samtools merge "${MERGED_BAM}" "${INDIV_DIR}"/*."${CELLTYPE}.${AGE}.${GENE}.bam"
    samtools index "${MERGED_BAM}"

    # Step 3: calculate subsampling fraction
    TOTAL_READS=$(samtools idxstats "${MERGED_BAM}" \
        | cut -f3 \
        | awk 'BEGIN{total=0} {total+=$1} END{print total}')
    FRACTION=$(echo "scale=6; ${SAMPLE_READS} / ${TOTAL_READS}" | bc -l)

    echo "  Total reads: ${TOTAL_READS}"
    echo "  Fraction   : ${FRACTION}"

    # Step 4: subsample or copy
    SAMPLED_BAM="${OUT_DIR}/mergeAll.${CELLTYPE}.${AGE}.${GENE}.sampled${SAMPLE_READS}.bam"

    if (( $(echo "${FRACTION} < 1" | bc -l) )); then
        samtools view -s "${FRACTION}" --subsample-seed 1 "${MERGED_BAM}" -o "${SAMPLED_BAM}"
    else
        cp "${MERGED_BAM}" "${SAMPLED_BAM}"
    fi

    samtools index "${SAMPLED_BAM}"

    # Step 5: write to TSV
    printf "%s\t%s\n" "${AGE}" "${SAMPLED_BAM}" >> "${TMP_TSV}"

    echo "  Done: ${AGE}"

done

mv "${TMP_TSV}" "${OUTPUT_TSV}"

echo ""
echo "========================================"
echo "  Finished! Output: ${OUTPUT_TSV}"
echo "========================================"


bash run_ggsashimi_age.sh "chr19:3430518-3488869" "NFIC" "CD8GZMH" 1000

ggsashimi.py \
    --bam         NFIC.tsv \
    --coordinates chr19:3455412-3463975 \
    --out-prefix  sashimi.NFIC \
    --gtf         gencode.v32.gtf \
    --palette     palette.age.txt



