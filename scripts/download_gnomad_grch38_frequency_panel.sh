#!/usr/bin/env bash
set -euo pipefail

# Build an ancestry panel from the public gnomAD v4.1 GRCh38 genome VCFs.
# gnomAD's genetic-ancestry groups are used as the five compatible columns:
# AFR, NFE (EUR proxy), EAS, SAS, and AMR. This is an ancestry-frequency
# panel, not a set of 1000 Genomes super-population frequencies.
#
# The source VCFs are queried through their public tabix indexes, so the full
# multi-gigabyte VCFs are not retained locally. Set GNOMAD_BASE_URL to the AWS
# or Azure mirror if GCS is not reachable from the cluster.

out_dir=${1:-data/ancestry}
min_af=${MIN_AF:-0.005}
max_af=${MAX_AF:-0.995}
min_spread=${MIN_SPREAD:-0.20}
min_an=${MIN_AN:-1000}
release=${GNOMAD_RELEASE:-4.1}
base_url=${GNOMAD_BASE_URL:-https://storage.googleapis.com/gcp-public-data--gnomad/release/${release}/vcf/genomes}
panel="$out_dir/gnomad_genomes_v${release}_grch38_pass_an${min_an}_ancestry.tsv"
shard_dir="$out_dir/gnomad_genomes_v${release}_grch38_pass_an${min_an}_shards"
mkdir -p "$shard_dir"

command -v bcftools >/dev/null || { echo 'bcftools is required' >&2; exit 127; }

for chromosome in $(seq 1 22); do
    shard="$shard_dir/chr${chromosome}.tsv"
    partial="$shard.part"
    if [[ -s "$shard" ]]; then
        printf 'Using completed chromosome %s shard: %s\n' "$chromosome" "$shard" >&2
        continue
    fi

    source="$base_url/gnomad.genomes.v${release}.sites.chr${chromosome}.vcf.bgz"
    rm -f "$shard" "$partial"

    # bcftools emits comma-separated values for Number=A fields. Restricting
    # REF and ALT to one base makes the first AF value unambiguous.
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AF_afr\t%INFO/AF_nfe\t%INFO/AF_eas\t%INFO/AF_sas\t%INFO/AF_amr\t%INFO/AN_afr\t%INFO/AN_nfe\t%INFO/AN_eas\t%INFO/AN_sas\t%INFO/AN_amr\n' "$source" \
          | awk -F '\t' -v OFS='\t' -v lo="$min_af" -v hi="$max_af" -v spread="$min_spread" -v minan="$min_an" '
          $3 ~ /^[ACGT]$/ && $4 ~ /^[ACGT]$/ && $3 != $4 && $5 == "PASS" &&
          $6 ~ /^[0-9.eE+-]+$/ && $7 ~ /^[0-9.eE+-]+$/ && $8 ~ /^[0-9.eE+-]+$/ &&
          $9 ~ /^[0-9.eE+-]+$/ && $10 ~ /^[0-9.eE+-]+$/ &&
          $11 ~ /^[0-9]+$/ && $12 ~ /^[0-9]+$/ && $13 ~ /^[0-9]+$/ &&
          $14 ~ /^[0-9]+$/ && $15 ~ /^[0-9]+$/ {
              mx=$6; mn=$6
              for (i=7; i<=10; i++) { if ($i > mx) mx=$i; if ($i < mn) mn=$i }
              if (mx >= lo && mn <= hi && mx-mn >= spread &&
                  $6 >= 0 && $6 <= 1 && $7 >= 0 && $7 <= 1 &&
                  $8 >= 0 && $8 <= 1 && $9 >= 0 && $9 <= 1 && $10 >= 0 && $10 <= 1 &&
                  $11 >= minan && $12 >= minan && $13 >= minan && $14 >= minan && $15 >= minan)
                  print $1,$2,$3,$4,$6,$7,$8,$9,$10
          }' > "$partial"
    [[ -s "$partial" ]] || { echo "Empty chromosome $chromosome shard" >&2; exit 1; }
    mv "$partial" "$shard"
done

printf 'chrom\tpos\tref\talt\tAFR\tNFE\tEAS\tSAS\tAMR\n' > "$panel.sorted"
for chromosome in $(seq 1 22); do
    cat "$shard_dir/chr${chromosome}.tsv"
done >> "$panel.sorted"
mv "$panel.sorted" "$panel"

printf 'Panel written: %s\n' "$panel"
printf 'Markers: '
tail -n +2 "$panel" | wc -l
printf 'Source: gnomAD genomes v%s, GRCh38 genetic-ancestry groups.\n' "$release"
printf 'Filters: PASS min_group_AN=%s min_af=%s max_af=%s min_population_spread=%s\n' "$min_an" "$min_af" "$max_af" "$min_spread"
