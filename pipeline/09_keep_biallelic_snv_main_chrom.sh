#!/bin/bash
# Keep biallelic SNVs, then keep chr1-22, X and Y.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate variant-calling

RNA_IN=${ROOT}/results/vcf
DNA_DIR=${ROOT}/singlecell_data
R5=${ROOT}/results/rna_step5_snv
D5=${ROOT}/results/dna_step5_snv
R6=${ROOT}/results/rna_step6_mainchr
D6=${ROOT}/results/dna_step6_mainchr
mkdir -p "$R5" "$D5" "$R6" "$D6"

SUM5=${ROOT}/results/step5_summary.tsv
SUM6=${ROOT}/results/step6_summary.tsv
printf "type\tcell\tinput\tafter_step5\n" > "$SUM5"
printf "type\tcell\tafter_step5\tafter_step6\n" > "$SUM6"

MAINCHR='^chr([1-9]|1[0-9]|2[0-2]|X|Y)$'

# Keep biallelic SNVs only
do_step5() {
  in="$1"; out="$2"; type="$3"; base="$4"
  # resume (added 2026-08-15): skip files already produced, so new cells can be added
  # without redoing the validated old ones.
  if [ -s "$out" ]; then echo "STEP5 skip (done) $type $base"; return; fi
  bcftools view -v snps -m2 -M2 "$in" -Oz -o "$out" 2>/dev/null
  tabix -p vcf "$out" 2>/dev/null
  tot=$(bcftools view -H "$in"  | wc -l)
  a5=$(bcftools view -H "$out" | wc -l)
  printf "%s\t%s\t%s\t%s\n" "$type" "$base" "$tot" "$a5" >> "$SUM5"
  echo "STEP5 $type $base: input=$tot  step5=$a5"
}

echo ">>> biallelic SNVs"
echo "RNA (35)"
find "$RNA_IN" -name "*.filtered.vcf.gz" | sort | while read f; do
  base=$(basename "$f" .filtered.vcf.gz)
  do_step5 "$f" "$R5/${base}.snv.vcf.gz" "RNA" "$base"
done
echo "DNA (35)"
find "$DNA_DIR" -name "*.vcf.gz" | sort | while read f; do
  base=$(basename "$f" .vcf.gz)
  do_step5 "$f" "$D5/${base}.snv.vcf.gz" "DNA" "$base"
done
echo ">>> step 5 done. Output: $R5 (RNA), $D5 (DNA). Summary: $SUM5"

# Keep the main chromosomes only
do_step6() {
  in="$1"; out="$2"; type="$3"; base="$4"
  # resume (added 2026-08-15) - see note in do_step5
  if [ -s "$out" ]; then echo "STEP6 skip (done) $type $base"; return; fi
  bcftools view -h "$in" 2>/dev/null > /tmp/hdr.$$
  bcftools view -H "$in" 2>/dev/null \
    | awk -v re="$MAINCHR" 'BEGIN{FS=OFS="\t"} $1 ~ re {print}' > /tmp/body.$$
  cat /tmp/hdr.$$ /tmp/body.$$ | bgzip > "$out"
  rm -f /tmp/hdr.$$ /tmp/body.$$
  tabix -p vcf "$out" 2>/dev/null
  a5=$(bcftools view -H "$in"  | wc -l)
  a6=$(bcftools view -H "$out" | wc -l)
  printf "%s\t%s\t%s\t%s\n" "$type" "$base" "$a5" "$a6" >> "$SUM6"
  echo "STEP6 $type $base: step5=$a5  step6=$a6"
}

echo ""
echo ">>> chr1-22, X, Y"
echo "RNA (35)"
find "$R5" -name "*.snv.vcf.gz" | sort | while read f; do
  base=$(basename "$f" .snv.vcf.gz)
  do_step6 "$f" "$R6/${base}.mainchr.vcf.gz" "RNA" "$base"
done
echo "DNA (35)"
find "$D5" -name "*.snv.vcf.gz" | sort | while read f; do
  base=$(basename "$f" .snv.vcf.gz)
  do_step6 "$f" "$D6/${base}.mainchr.vcf.gz" "DNA" "$base"
done
echo ">>> step 6 done. Output: $R6 (RNA), $D6 (DNA). Summary: $SUM6"

echo ""
echo ">>> all done (steps 5 + 6)."
