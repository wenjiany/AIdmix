#!/usr/bin/env bash
set -euo pipefail

# Audit only files already tracked by this clean public repository. This is a
# guardrail, not a substitute for human review before changing visibility.
forbidden_path='(^|/)(data|results|logs|work|tmp|hippafy_xref)(/|$)|[.](bam|bai|cram|crai|sam|vcf|bcf|tbi|csi|fa|fasta|fna|fai|gzi|dict|bed|xlsx|xls|csv|tsv|list|log|out|err)$'
if git ls-files | grep -Ei "$forbidden_path"; then
    printf 'ERROR: tracked data-like or result-like path detected\n' >&2
    exit 1
fi

forbidden_content='/research/|/home/|authorized_apps|rgs01|yanggrp|TTCRC|ONTDNA|SUBJECT-ID|hippafy_xref'
if git grep -InE "$forbidden_content" -- +    ':!.gitignore' +    ':!scripts/audit_public_tree.sh'; then
    printf 'ERROR: tracked content contains a local or identifying marker\n' >&2
    exit 1
fi

large_files=$(git ls-files -z | xargs -0 -r stat -c '%s %n' | awk '$1 > 1000000')
if [[ -n "$large_files" ]]; then
    printf 'ERROR: tracked files larger than 1 MB:\n%s\n' "$large_files" >&2
    exit 1
fi

printf 'Public-tree audit passed: %s tracked files\n' "$(git ls-files | wc -l)"
