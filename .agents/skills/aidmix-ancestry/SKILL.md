---
name: aidmix-ancestry
description: Prepare validated ancestry-frequency panels and estimate genetic ancestry from indexed BAM or CRAM alignments with AIdmix, including hg38/T2T liftover, pilot and cohort runs, aggregation, comparison, and QC. Use for AIdmix workflows or methodological comparisons with iAdmix.
license: MIT
metadata:
  author: "Wenjian Yang"
  version: "0.1.0"
  repository: "https://github.com/wenjiany/AIdmix"
---

# AIdmix ancestry estimation

AIdmix is an AI-assisted implementation of a conventional genotype-likelihood
mixture model. The inference algorithm itself is not an AI model. Acknowledge
iAdmix as methodological inspiration when documenting or publishing the work.

The workflow requires Python 3.10+, `samtools`, an indexed BAM or CRAM, its
matching reference FASTA, and a compatible ancestry-frequency panel.

## Locate the implementation

Prefer the installed console commands and resolve them with `command -v`:

- `aidmix`
- `aidmix-prepare-panel`
- `aidmix-combine`

If they are unavailable, use `AIDMIX_ROOT` when defined or locate a checkout of
[wenjiany/AIdmix](https://github.com/wenjiany/AIdmix). Do not download or
install software without the user's authorization. Keep the version-controlled
AIdmix repository authoritative rather than copying implementation code into
the skill.

Read [references/methodology.md](references/methodology.md) when explaining,
reviewing, or changing the statistical model, gradient, Hessian, optimizer, or
read-error assumptions.

## Non-negotiable input checks

The alignment, panel, and reference FASTA must describe the same assembly and
compatible contig names. Matching chromosome lengths is insufficient: use the
exact FASTA associated with the BAM/CRAM whenever possible.

Every normalized panel row must be reconciled against that FASTA:

- keep the row when panel REF equals FASTA REF;
- when panel ALT equals FASTA REF, swap REF/ALT and complement every population
  ALT frequency as `1 - frequency`;
- exclude and report rows for which neither allele matches;
- fail the workflow when marker retention or observed-site counts are
  unexpectedly low.

The current estimator performs this reconciliation at load time. Preserve it
when adapting the code. For a newly lifted panel, independently verify that REF
mismatches are zero before a cohort run.

## Panel formats and interpretation

Normalized format:

```text
chrom  pos  ref  alt  population_1  population_2 ...
```

Native iAdmix format is also accepted by the estimator:

```text
#chrom position rsid A1 A2 population_1 population_2 ...
```

Native iAdmix frequencies describe A1. Orient A1/A2 against the source FASTA
before liftover; do not merely rename A1 and A2 to REF and ALT.

Common panels are not interchangeable:

- Omni: five populations (`european`, `native`, `eastasian`, `southasian`,
  `african`).
- SNP6: six populations, splitting European into `north.european` and
  `south.european`.
- gnomAD: ancestry groups such as AFR, EUR/NFE, EAS, SAS, and AMR. AMR is not
  equivalent to the iAdmix Native component.

Record panel source, version, build, population columns, filtering, liftover
chains, target FASTA, marker count, and REF-validation count with every run.

## Prepare or lift a panel

`aidmix-prepare-panel` expects normalized input. With hg38 input and T2T
output, require unique reciprocal liftover and target-reference validation:

```bash
aidmix-prepare-panel \
  --input panel.hg38.normalized.tsv \
  --output panel.t2t.tsv \
  --reference /path/to/T2T.fa \
  --chain /path/to/hg38ToT2T.chain.gz \
  --reverse-chain /path/to/T2TToHg38.chain.gz \
  --min-spread 0.20 \
  --min-distance 500
```

Validate source REF alleles against the source FASTA before liftover and output
REF alleles against the target FASTA afterward. Keep superseded invalid panels
clearly isolated so they cannot be selected accidentally.

## Run a pilot

Measure one representative alignment before a cohort. Example:

```bash
aidmix \
  --panel /path/to/panel.tsv \
  --reference /path/to/matching.fa \
  --alignment /path/to/sample.cram \
  --output /path/to/sample.ancestry.tsv \
  --min-mapq 20 \
  --min-baseq 10 \
  --min-sites 100 \
  --min-chromosomes 10 \
  --bootstrap 0
```

Typical hg38 RNA-seq/iAdmix-compatible filtering uses MAPQ 30 and base quality
13; ONT WGS pilots have used MAPQ 20 and base quality 10. Treat these as
starting points, not universal platform requirements. The default
`--max-depth 10000` is a high safety ceiling; lowering it is not assumed to
improve runtime and should be benchmarked on the actual alignment workload.

Before scaling, require:

- nonzero, plausible observed-site and read counts;
- broad autosomal coverage for WGS;
- finite likelihoods and successful KKT convergence;
- ancestry proportions summing to approximately 100%;
- empirically measured runtime and memory.

## Run a cohort

Use the local batch scheduler conventions when available. Size chunks from the
pilot so jobs are long enough to schedule efficiently. Each sample must write
one atomic result, and chunk runners must skip complete outputs, continue after
a sample failure, print explicit failure markers, and print a completion marker.

Never edit a shared runner while submitted jobs are pending or running.
Monitor output files rather than only scheduler state. Before declaring
completion, reconcile every input alignment against its expected output and
inspect scheduler termination/memory messages in addition to explicit failures.

Keep bootstrap disabled for the initial point-estimate cohort. Add
chromosome-block bootstrap only after point-estimate QC and resource sizing.

## Combine and compare

```bash
aidmix-combine \
  --input-dir /path/to/per_sample \
  --output /path/to/cohort.ancestry.tsv
```

Compare panels only after harmonizing allele orientation and population
definitions. Report dominant-component agreement, per-component correlation,
absolute percentage-point differences, likelihood/convergence status, observed
sites, and chromosome coverage. Near-zero components can have weak correlation
despite negligible absolute differences.

## Privacy and sharing

Do not commit or package alignments, sample lists, patient identifiers, linkage
tables, cohort outputs, scheduler logs containing identifiers, or institutional
paths. The public repository and this skill should contain only source code,
portable instructions, tests, and non-identifying methodology.
