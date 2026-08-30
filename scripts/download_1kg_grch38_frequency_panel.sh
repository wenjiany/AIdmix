#!/usr/bin/env bash
set -euo pipefail

# Stream official IGSR/EBI phase-3 GRCh38 VCFs and retain only the compact
# five-superpopulation frequency panel used by AIdmix.
# No multi-gigabyte source VCF is retained locally.

out_dir=${1:-data/ancestry}
min_af=${MIN_AF:-0.005}
max_af=${MAX_AF:-0.995}
min_spread=${MIN_SPREAD:-0.20}
mkdir -p "$out_dir"
base='https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/phase3_liftover_nygc_dir'
manifest="$out_dir/phase3_grch38_manifest.tsv"
panel="$out_dir/1kg_phase3_grch38_superpop.tsv"
shard_dir="$out_dir/1kg_phase3_grch38_shards"
mkdir -p "$shard_dir"

curl --fail --silent --show-error --location --retry 6 --retry-all-errors --retry-delay 15 \
    "$base/phase3.crossmap.GRCh38.07302021.manifest.tsv" -o "$manifest"
printf 'chrom\tpos\tref\talt\tAFR\tEUR\tEAS\tSAS\tAMR\n' > "$panel"

for chromosome in $(seq 1 22); do
    shard="$shard_dir/chr${chromosome}.tsv"
    if [[ -s "$shard" ]]; then
        printf 'Using completed chromosome %s shard: %s\n' "$chromosome" "$shard" >&2
        continue
    fi

    source="$base/phase3.chr${chromosome}.GRCh38.GT.crossmap.vcf.gz"
    temp_root="${TMPDIR:-/tmp}/aidmix_1kg_chr${chromosome}"
    temp_vcf="${temp_root}.vcf.gz.part"
    temp_tbi="${temp_vcf}.tbi.part"
    rm -f "$shard"
    curl --fail --silent --show-error --location --retry 6 --retry-all-errors --retry-delay 15 \
        --continue-at - "$source" -o "$temp_vcf"
    curl --fail --silent --show-error --location --retry 6 --retry-all-errors --retry-delay 15 \
        --continue-at - "$source.tbi" -o "$temp_tbi"
    # Validate against the official manifest. Some files have inconsistent
    # VCF contig metadata (for example chr2 records with an ID=2 header), so
    # bcftools is not used for extraction below.
    expected_md5=$(awk -v f="$(basename "$source")" '$1 == f {print $3}' "$manifest")
    actual_md5=$(md5sum "$temp_vcf" | awk '{print $1}')
    if [[ -z "$expected_md5" || "$actual_md5" != "$expected_md5" ]] || ! gzip -t "$temp_vcf" >/dev/null 2>&1; then
        printf 'Invalid or incomplete transfer for chromosome %s; restarting\n' "$chromosome" >&2
        rm -f "$temp_vcf" "$temp_tbi"
        curl --fail --silent --show-error --location --retry 6 --retry-all-errors --retry-delay 15 \
            "$source" -o "$temp_vcf"
        curl --fail --silent --show-error --location --retry 6 --retry-all-errors --retry-delay 15 \
            "$source.tbi" -o "$temp_tbi"
        actual_md5=$(md5sum "$temp_vcf" | awk '{print $1}')
        [[ "$actual_md5" == "$expected_md5" ]] || { printf 'MD5 mismatch for chromosome %s\n' "$chromosome" >&2; exit 1; }
        gzip -t "$temp_vcf"
    fi
    gzip -cd "$temp_vcf" \
      | awk -F '\t' -v OFS='\t' -v lo="$min_af" -v hi="$max_af" -v spread="$min_spread" '
          /^#/ {next}
          length($4) != 1 || length($5) != 1 || $4 !~ /^[ACGT]$/ || $5 !~ /^[ACGT]$/ || $5 ~ /,/ {next}
          {
              afr=eur=eas=sas=amr=""; n=split($8, info, ";")
              for (i=1; i<=n; i++) {
                  split(info[i], kv, "=")
                  if (kv[1] == "AFR_AF") afr=kv[2]
                  else if (kv[1] == "EUR_AF") eur=kv[2]
                  else if (kv[1] == "EAS_AF") eas=kv[2]
                  else if (kv[1] == "SAS_AF") sas=kv[2]
                  else if (kv[1] == "AMR_AF") amr=kv[2]
              }
              if (afr == "" || eur == "" || eas == "" || sas == "" || amr == "") next
              mx=afr; mn=afr; v=eur; if (v > mx) mx=v; if (v < mn) mn=v
              v=eas; if (v > mx) mx=v; if (v < mn) mn=v
              v=sas; if (v > mx) mx=v; if (v < mn) mn=v
              v=amr; if (v > mx) mx=v; if (v < mn) mn=v
              if (mx >= lo && mn <= hi && mx-mn >= spread && afr >= 0 && afr <= 1 && eur >= 0 && eur <= 1 && eas >= 0 && eas <= 1 && sas >= 0 && sas <= 1 && amr >= 0 && amr <= 1) print $1,$2,$4,$5,afr,eur,eas,sas,amr
          }' > "$shard"
    if [[ ! -s "$shard" ]]; then
        printf 'Empty chromosome %s shard; preserving source for diagnosis\n' "$chromosome" >&2
        exit 1
    fi
    rm -f "$temp_vcf" "$temp_tbi"
done

# Assemble completed per-chromosome shards, remove duplicate coordinates,
# retain the first ALT projection, and sort. A missing shard is an error rather
# than silently producing an incomplete panel.
for chromosome in $(seq 1 22); do
    test -s "$shard_dir/chr${chromosome}.tsv" || {
        printf 'Missing chromosome %s shard; panel not assembled\n' "$chromosome" >&2
        exit 1
    }
done
{
    printf 'chrom\tpos\tref\talt\tAFR\tEUR\tEAS\tSAS\tAMR\n'
    cat "$shard_dir"/chr{1..22}.tsv
} | awk 'NR == 1 {print; next} !seen[$1 FS $2]++' \
  | sort -k1,1V -k2,2n > "$panel.sorted"
mv "$panel.sorted" "$panel"

printf 'Panel written: %s\nMarkers: ' "$panel"
tail -n +2 "$panel" | wc -l
printf 'Filters: min_af=%s max_af=%s min_population_spread=%s\n' "$min_af" "$max_af" "$min_spread"
