# AIdmix

AIdmix estimates genetic-ancestry mixture proportions directly from
coordinate-sorted, indexed BAM or CRAM alignments. It combines read-level
genotype likelihoods with population allele frequencies and maximizes the
joint likelihood over ancestry proportions. It does not require hard genotype
calls and is not tied to a sequencing platform or reference assembly.

The inference algorithm is a conventional statistical likelihood model, not
an AI model. AIdmix is similar in overall statistical purpose and was
methodologically inspired by the [iAdmix project](https://github.com/eliorav/iAdmix),
which estimates admixture coefficients from genotype or sequence data. AIdmix
is an independent implementation and is not affiliated with iAdmix.

The two projects should not be assumed to produce identical estimates: panel
definitions, read/genotype-likelihood construction, reference assemblies, and
optimization details may differ.

## Requirements

- Python 3.10 or newer;
- `samtools` available on `PATH`;
- a coordinate-sorted and indexed BAM or CRAM;
- the exact matching reference FASTA and index; and
- an ancestry-frequency panel on the same assembly.

Install the tagged release directly from GitHub:

```bash
python -m pip install 'aidmix-ancestry @ git+https://github.com/wenjiany/AIdmix.git@v0.1.0'
```

Or clone the repository and install it for development:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[test]'
```

## Panel format

The normalized tab-separated format is:

```text
chrom  pos  ref  alt  population_1  population_2 ...
```

Positions are one-based and every population column contains the ALT-allele
frequency. Native iAdmix panels are also accepted, but must be oriented against
the supplied reference FASTA at runtime.

For normalized panels, AIdmix verifies every REF allele against the FASTA. If
the FASTA carries the panel ALT as REF, AIdmix swaps REF/ALT and complements
all frequencies as `1-p`. Markers for which neither allele matches are
excluded and reported.

## Prepare or lift a panel

```bash
aidmix-prepare-panel \
  --input panel.source.tsv \
  --output panel.target.tsv \
  --reference /path/to/target.fa \
  --chain /path/to/sourceToTarget.chain.gz \
  --reverse-chain /path/to/targetToSource.chain.gz \
  --min-spread 0.20 \
  --min-distance 500
```

Omit both chain arguments when the input already uses the target assembly.
The scripts in `scripts/` construct possible GRCh38 frequency sources from
1000 Genomes or gnomAD. Generated panels are intentionally excluded from Git.
Population labels and definitions from different sources are not necessarily
interchangeable.

## Estimate ancestry

```bash
aidmix \
  --panel /path/to/panel.tsv \
  --reference /path/to/reference.fa \
  --alignment /path/to/sample.bam \
  --output /path/to/sample.ancestry.tsv \
  --min-mapq 20 \
  --min-baseq 10 \
  --min-sites 100 \
  --min-chromosomes 10 \
  --bootstrap 0
```

Multiple alignments can be supplied by repeating `--alignment` or by using
`--alignment-list`. When a target BED is supplied, AIdmix reports background,
target, and downweighted joint estimates. Without a target BED, observed sites
are treated as background.

## Combine one-sample outputs

```bash
aidmix-combine \
  --input-dir /path/to/per-sample-results \
  --output /path/to/cohort.ancestry.tsv
```

The combiner requires exactly one row per input and identical headers.

## Method and interpretation

See [the methodology document](docs/methodology.md) for the read likelihood,
genotype marginalization, ancestry likelihood, analytic gradient and Hessian,
optimization, and KKT validation.

Ancestry estimates are coordinates relative to the selected reference panel.
They are not race or ethnicity, and labels from different panels should not be
treated as equivalent without validation.

## Testing

```bash
python -m unittest discover -s tests
```

## Privacy

This repository contains source code, documentation, and synthetic tests only.
Do not commit alignments, reference genomes, frequency panels, sample lists,
cohort outputs, linkage tables, scheduler logs, identifiers, or institutional
filesystem paths. Run `scripts/audit_public_tree.sh` before publishing.
