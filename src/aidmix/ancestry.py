#!/usr/bin/env python3
"""Estimate ancestry proportions from coordinate-sorted BAM or CRAM files.

The input panel is a tab-separated marker table with columns:
chrom, pos, ref, alt, followed by two or more population-frequency columns.
Positions are 1-based and frequencies refer to the ALT allele.
Reads are evaluated directly and diploid genotypes are marginalized without
requiring hard genotype calls.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import subprocess
import sys
from pathlib import Path
from dataclasses import dataclass

import numpy as np
import pysam
from scipy.optimize import minimize

BASES = frozenset("ACGT")
AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
BASE_COLUMNS = ("chrom", "pos", "ref", "alt")


@dataclass(frozen=True)
class Marker:
    chrom: str
    pos: int
    ref: str
    alt: str
    freq: tuple[float, ...]


@dataclass
class SiteObs:
    marker: Marker
    region: str
    ll: tuple[float, float, float]
    reads: int


def canonical_chrom(value: str) -> str:
    return value if value.startswith("chr") else "chr" + value


def read_bed(path: str | None, buffer: int) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = {}
    if not path:
        return out
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) < 3:
                continue
            chrom, start, end = canonical_chrom(f[0]), int(f[1]), int(f[2])
            out.setdefault(chrom, []).append((max(0, start - buffer), end + buffer))
    for chrom in out:
        out[chrom].sort()
    return out


def in_bed(bed: dict[str, list[tuple[int, int]]], chrom: str, pos0: int) -> bool:
    return any(start <= pos0 < end for start, end in bed.get(chrom, ()))


def read_panel(path: str, min_af: float, max_af: float,
               reference: str | None = None) -> tuple[list[Marker], tuple[str, ...]]:
    """Read normalized panels or native iAdmix ``#chrom position ...`` panels."""
    markers: list[Marker] = []
    reoriented = 0
    reference_mismatches = 0
    with open(path) as handle:
        header_line = next((line for line in handle if line.strip()), "")
        header = header_line.lstrip("#").split()
        legacy = {"chrom", "position", "rsid", "A1", "A2"}.issubset(header)
        normalized = set(BASE_COLUMNS).issubset(header)
        if not legacy and not normalized:
            raise SystemExit("Panel must use chrom/pos/ref/alt format or native iAdmix format")
        if legacy and reference is None:
            raise SystemExit("--reference is required to orient a native iAdmix panel")
        if legacy:
            populations = tuple(header[header.index("A2") + 1:])
        else:
            populations = tuple(x for x in header if x not in BASE_COLUMNS)
        if len(populations) < 2:
            raise SystemExit("Panel must contain at least two population-frequency columns")
        indices = {name: header.index(name) for name in header}
        fasta = pysam.FastaFile(reference) if reference else None
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            try:
                chrom = canonical_chrom(fields[indices["chrom"]])
                pos = int(fields[indices["position"] if legacy else indices["pos"]])
                if legacy:
                    allele1 = fields[indices["A1"]].upper()
                    allele2 = fields[indices["A2"]].upper()
                    freq = tuple(float(fields[indices[p]]) for p in populations)
                    ref = fasta.fetch(chrom, pos - 1, pos).upper()
                    if ref == allele1:
                        alt = allele2
                        freq = tuple(1.0 - value for value in freq)
                    elif ref == allele2:
                        alt = allele1
                    else:
                        continue
                else:
                    ref = fields[indices["ref"]].upper()
                    alt = fields[indices["alt"]].upper()
                    freq = tuple(float(fields[indices[p]]) for p in populations)
                    if fasta is not None:
                        actual_ref = fasta.fetch(chrom, pos - 1, pos).upper()
                        if actual_ref == alt:
                            ref, alt = alt, ref
                            freq = tuple(1.0 - value for value in freq)
                            reoriented += 1
                        elif actual_ref != ref:
                            reference_mismatches += 1
                            continue
            except (TypeError, ValueError):
                continue
            if chrom not in AUTOSOMES or len(ref) != 1 or len(alt) != 1 or ref not in BASES or alt not in BASES or ref == alt:
                continue
            if any(not 0 <= q <= 1 for q in freq) or max(freq) < min_af or min(freq) > max_af:
                continue
            markers.append(Marker(chrom, pos, ref, alt, freq))
    if reoriented or reference_mismatches:
        print(
            f"Panel/reference reconciliation: reoriented={reoriented} "
            f"unresolved={reference_mismatches}",
            file=sys.stderr,
        )
    unique = {(m.chrom, m.pos): m for m in markers}
    return sorted(unique.values(), key=lambda m: (int(m.chrom[3:]), m.pos)), populations


def logsumexp(values: tuple[float, ...]) -> float:
    peak = max(values)
    return peak + math.log(sum(math.exp(x - peak) for x in values))


def base_error(q: int) -> float:
    return min(0.25, max(0.005, 10.0 ** (-max(1, q) / 10.0)))


def read_log_likelihood(ref: str, alt: str, observed: str, q: int, genotype: int) -> float:
    e = base_error(q)
    alt_fraction = genotype / 2.0
    if observed == ref:
        probability = (1 - alt_fraction) * (1 - e) + alt_fraction * e / 3
    elif observed == alt:
        probability = alt_fraction * (1 - e) + (1 - alt_fraction) * e / 3
    else:
        probability = e / 3
    return math.log(max(probability, 1e-12))


def pileup_sites(cram: str, reference: str, markers: list[Marker], target_bed: dict[str, list[tuple[int, int]]],
                 min_mapq: int, min_baseq: int, max_depth: int) -> list[SiteObs]:
    by_coord = {(m.chrom, m.pos): m for m in markers}
    sites_bed = Path(os.environ.get("TMPDIR", "/tmp")) / f"aim_sites_{os.getpid()}.bed"
    try:
        with sites_bed.open("w") as out:
            for m in markers:
                out.write(f"{m.chrom}\t{m.pos - 1}\t{m.pos}\n")
        cmd = ["samtools", "mpileup", "-B", "-f", reference, "-l", str(sites_bed),
               "-q", str(min_mapq), "-Q", str(min_baseq), "-d", str(max_depth), cram]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        observations: list[SiteObs] = []
        assert proc.stdout is not None
        for line in proc.stdout:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            marker = by_coord.get((canonical_chrom(f[0]), int(f[1])))
            if marker is None:
                continue
            region = "target" if in_bed(target_bed, marker.chrom, marker.pos - 1) else "background"
            bases, qualities = f[4], f[5]
            ll = [0.0, 0.0, 0.0]
            n = 0
            i = 0
            qi = 0
            while i < len(bases) and qi < len(qualities):
                c = bases[i]
                if c == "^":
                    i += 2
                    continue
                if c == "$":
                    i += 1
                    continue
                if c in "+-":
                    i += 1
                    j = i
                    while j < len(bases) and bases[j].isdigit():
                        j += 1
                    if j > i:
                        i = j + int(bases[i:j])
                    continue
                q = ord(qualities[qi]) - 33
                qi += 1
                if c in ".,":
                    observed = marker.ref
                elif c.upper() in BASES:
                    observed = c.upper()
                else:
                    i += 1
                    continue
                if observed not in (marker.ref, marker.alt):
                    i += 1
                    continue
                for g in range(3):
                    ll[g] += read_log_likelihood(marker.ref, marker.alt, observed, q, g)
                n += 1
                i += 1
            if n:
                observations.append(SiteObs(marker, region, tuple(ll), n))
        stderr = proc.stderr.read() if proc.stderr is not None else ""
        rc = proc.wait()
        if rc:
            raise RuntimeError(f"samtools mpileup failed for {cram}: {stderr[-500:]}")
        return observations
    finally:
        sites_bed.unlink(missing_ok=True)


def softmax(x: np.ndarray) -> np.ndarray:
    z = np.r_[x, 0.0]
    z -= np.max(z)
    e = np.exp(z)
    return e / e.sum()


def score_data(observations: list[SiteObs], target_weight: float,
               region: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    usable = [o for o in observations if region == "joint" or o.region == region]
    if not usable:
        n_pops = len(observations[0].marker.freq) if observations else 0
        return np.empty((0, n_pops)), np.empty((0, 3)), np.empty(0)
    frequencies = np.asarray([o.marker.freq for o in usable], dtype=float)
    log_likelihoods = np.asarray([o.ll for o in usable], dtype=float)
    weights = np.asarray([
        target_weight if region == "joint" and o.region == "target" else 1.0
        for o in usable
    ])
    return frequencies, log_likelihoods, weights


def score_ancestry_arrays(ancestry: np.ndarray, frequencies: np.ndarray,
                          log_likelihoods: np.ndarray,
                          weights: np.ndarray) -> tuple[float, np.ndarray]:
    """Return vectorized log likelihood and gradient over ancestry proportions."""
    raw_q = frequencies @ ancestry
    q = np.clip(raw_q, 1e-6, 1 - 1e-6)
    priors = np.column_stack(((1 - q) ** 2, 2 * q * (1 - q), q ** 2))
    prior_derivatives = np.column_stack((-2 * (1 - q), 2 - 4 * q, 2 * q))
    peaks = np.max(log_likelihoods, axis=1)
    scaled_likelihoods = np.exp(log_likelihoods - peaks[:, None])
    marginal = np.sum(scaled_likelihoods * priors, axis=1)
    marginal = np.maximum(marginal, 1e-300)
    value = float(np.sum(weights * (peaks + np.log(marginal))))

    dlog_dq = np.sum(scaled_likelihoods * prior_derivatives, axis=1) / marginal
    clipped = (raw_q <= 1e-6) | (raw_q >= 1 - 1e-6)
    dlog_dq[clipped] = 0.0
    ancestry_gradient = frequencies.T @ (weights * dlog_dq)
    return value, ancestry_gradient


def score_arrays(x: np.ndarray, frequencies: np.ndarray, log_likelihoods: np.ndarray,
                 weights: np.ndarray) -> tuple[float, np.ndarray]:
    """Return vectorized log likelihood and gradient over softmax logits."""
    ancestry = softmax(x)
    value, ancestry_gradient = score_ancestry_arrays(
        ancestry, frequencies, log_likelihoods, weights
    )
    centered = ancestry_gradient - float(np.dot(ancestry, ancestry_gradient))
    logit_gradient = ancestry[:-1] * centered[:-1]
    return value, logit_gradient


def score(x: np.ndarray, observations: list[SiteObs], target_weight: float, region: str) -> float:
    frequencies, log_likelihoods, weights = score_data(observations, target_weight, region)
    if not len(frequencies):
        return 0.0
    return score_arrays(x, frequencies, log_likelihoods, weights)[0]


def choose_fit_result(results: list[object]) -> object | None:
    finite = [result for result in results if math.isfinite(float(result.fun))]
    if not finite:
        return None
    best = min(finite, key=lambda result: result.fun)
    converged = [result for result in finite if result.success]
    if converged:
        best_converged = min(converged, key=lambda result: result.fun)
        tolerance = max(1e-6, 1e-8 * abs(float(best.fun)))
        if best_converged.fun <= best.fun + tolerance:
            return best_converged
    return best


def kkt_violation(ancestry: np.ndarray, gradient: np.ndarray,
                  n_sites: int, boundary_tolerance: float = 1e-6) -> float:
    """Return a scale-normalized KKT violation for simplex maximization."""
    scale = max(1.0, float(n_sites))
    scaled = gradient / scale
    active = ancestry > boundary_tolerance
    if not np.any(active):
        return float("inf")
    multiplier = float(np.mean(scaled[active]))
    active_violation = float(np.max(np.abs(scaled[active] - multiplier)))
    inactive_violation = float(np.max(np.maximum(scaled[~active] - multiplier, 0.0))) if np.any(~active) else 0.0
    return max(active_violation, inactive_violation)


def fit(observations: list[SiteObs], target_weight: float, region: str,
        n_pops: int | None = None) -> tuple[np.ndarray, float, bool]:
    n_pops = n_pops if n_pops is not None else (len(observations[0].marker.freq) if observations else 0)
    frequencies, log_likelihoods, weights = score_data(observations, target_weight, region)
    if not len(frequencies):
        return np.full(n_pops, np.nan), float("nan"), False
    candidates: list[tuple[float, np.ndarray, bool, float]] = []

    def add_candidate(ancestry: np.ndarray, value: float, success: bool) -> None:
        ancestry = np.clip(np.asarray(ancestry, dtype=float), 0.0, 1.0)
        ancestry /= ancestry.sum()
        value, gradient = score_ancestry_arrays(
            ancestry, frequencies, log_likelihoods, weights
        )
        violation = kkt_violation(ancestry, gradient, len(frequencies))
        candidates.append((value, ancestry, success, violation))

    # Logit-space optimization avoids explicit simplex constraints and handles
    # boundary optima without SLSQP's common false failure status.
    logit_starts = [np.zeros(n_pops - 1)]
    for index in range(n_pops):
        start = np.full(n_pops - 1, -12.0)
        if index < n_pops - 1:
            start[index] = 12.0
        logit_starts.append(start)
    for start in logit_starts:
        def logit_objective(x: np.ndarray) -> tuple[float, np.ndarray]:
            value, gradient = score_arrays(x, frequencies, log_likelihoods, weights)
            return -value, -gradient

        result = minimize(
            logit_objective,
            start,
            jac=True,
            method="L-BFGS-B",
            bounds=[(-30.0, 30.0)] * (n_pops - 1),
            options={"maxiter": 1000, "ftol": 1e-12, "gtol": 1e-8, "maxls": 50},
        )
        add_candidate(softmax(result.x), -float(result.fun), bool(result.success))

    # Keep constrained optimization as a fallback for unusual likelihood
    # surfaces or numerical issues in the logit parameterization.
    starts = [np.full(n_pops, 1 / n_pops)]
    for index in range(n_pops):
        start = np.zeros(n_pops)
        start[index] = 1.0
        starts.append(start)
    for start in starts:
        def objective(ancestry: np.ndarray) -> tuple[float, np.ndarray]:
            value, gradient = score_ancestry_arrays(
                ancestry, frequencies, log_likelihoods, weights
            )
            return -value, -gradient

        result = minimize(
            objective,
            start,
            jac=True,
            method="SLSQP",
            bounds=[(0.0, 1.0)] * n_pops,
            constraints={
                "type": "eq",
                "fun": lambda ancestry: float(np.sum(ancestry) - 1.0),
                "jac": lambda ancestry: np.ones_like(ancestry),
            },
            options={"maxiter": 500, "ftol": 1e-8},
        )
        add_candidate(np.asarray(result.x), -float(result.fun), bool(result.success))

    finite = [candidate for candidate in candidates if math.isfinite(candidate[0])]
    if not finite:
        return np.full(n_pops, np.nan), float("nan"), False
    valid = [candidate for candidate in finite if candidate[3] <= 1e-5]
    pool = valid if valid else finite
    value, ancestry, success, violation = max(pool, key=lambda candidate: candidate[0])
    converged = violation <= 1e-5
    return ancestry, value, converged


def bootstrap(observations: list[SiteObs], target_weight: float, region: str, reps: int, seed: int) -> np.ndarray:
    usable = [o for o in observations if region == "joint" or o.region == region]
    n_pops = len(observations[0].marker.freq) if observations else 0
    if len(usable) < 2 or reps <= 0:
        return np.empty((0, n_pops))
    rng = np.random.default_rng(seed)
    by_chrom = {}
    for obs in usable:
        by_chrom.setdefault(obs.marker.chrom, []).append(obs)
    chromosomes = sorted(by_chrom)
    results = []
    for _ in range(reps):
        selected = rng.choice(chromosomes, size=len(chromosomes), replace=True)
        sample = [o for chromosome in selected for o in by_chrom[chromosome]]
        results.append(fit(sample, target_weight, region)[0])
    return np.asarray(results)


def ancestry_one(cram: str, args: argparse.Namespace, markers: list[Marker], populations: tuple[str, ...],
                 bed: dict[str, list[tuple[int, int]]]) -> dict[str, object]:
    observations = pileup_sites(cram, args.reference, markers, bed, args.min_mapq, args.min_baseq, args.max_depth)
    alignment_name = Path(cram).name
    for suffix in (".T2T.haplotagged.cram", ".cram", ".bam"):
        if alignment_name.endswith(suffix):
            alignment_name = alignment_name[:-len(suffix)]
            break
    row: dict[str, object] = {"sample": alignment_name,
                              "observed_sites": len(observations), "observed_reads": sum(o.reads for o in observations),
                              "target_sites": sum(o.region == "target" for o in observations),
                              "background_sites": sum(o.region == "background" for o in observations),
                              "target_chromosomes": len({o.marker.chrom for o in observations if o.region == "target"}),
                              "background_chromosomes": len({o.marker.chrom for o in observations if o.region == "background"})}
    convergence = {}
    for region in ("background", "target", "joint"):
        ancestry, likelihood, converged = fit(
            observations, args.target_weight, region, n_pops=len(populations)
        )
        convergence[region] = converged
        boot = bootstrap(observations, args.target_weight, region, args.bootstrap, args.seed)
        row[f"{region}_loglik"] = f"{likelihood:.6f}" if math.isfinite(likelihood) else "."
        row[f"{region}_converged"] = "YES" if converged else "NO"
        for index, pop in enumerate(populations):
            row[f"{region}_{pop}_pct"] = f"{100 * ancestry[index]:.4f}" if math.isfinite(ancestry[index]) else "."
            row[f"{region}_{pop}_lo95"] = f"{100 * np.quantile(boot[:, index], .025):.4f}" if len(boot) else "."
            row[f"{region}_{pop}_hi95"] = f"{100 * np.quantile(boot[:, index], .975):.4f}" if len(boot) else "."
    if row["background_sites"] < args.min_sites or row["background_chromosomes"] < args.min_chromosomes:
        row["qc_status"] = "LOW_BACKGROUND_SIGNAL"
    elif not convergence["background"]:
        row["qc_status"] = "OPTIMIZER_WARNING"
    else:
        row["qc_status"] = "PASS"
    return row


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--panel", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--alignment", action="append")
    p.add_argument("--alignment-list")
    p.add_argument("--target-bed")
    p.add_argument("--target-buffer", type=int, default=1000)
    p.add_argument("--output", required=True)
    p.add_argument("--min-af", type=float, default=0.005)
    p.add_argument("--max-af", type=float, default=0.995)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--min-baseq", type=int, default=10)
    p.add_argument("--max-depth", type=int, default=10000)
    p.add_argument("--target-weight", type=float, default=0.25)
    p.add_argument("--bootstrap", type=int, default=100)
    p.add_argument("--seed", type=int, default=20260828)
    p.add_argument("--min-sites", type=int, default=100)
    p.add_argument("--min-chromosomes", type=int, default=10)
    args = p.parse_args()
    if not args.alignment and not args.alignment_list:
        p.error("provide --alignment or --alignment-list")
    alignments = list(args.alignment or [])
    if args.alignment_list:
        with open(args.alignment_list) as handle:
            alignments.extend(line.strip() for line in handle if line.strip() and not line.startswith("#"))
    markers, populations = read_panel(args.panel, args.min_af, args.max_af, args.reference)
    if not markers:
        raise SystemExit("No usable markers remain after panel filtering")
    bed = read_bed(args.target_bed, args.target_buffer)
    rows = [ancestry_one(alignment, args, markers, populations, bed) for alignment in alignments]
    fields = list(rows[0]) if rows else ["sample"]
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    print(f"Processed {len(rows)} samples at {len(markers)} panel markers; wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
