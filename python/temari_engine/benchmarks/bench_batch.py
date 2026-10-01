"""Reproducible scalar-vs-batch benchmark for the shipped F and factors tables.

Run from the repository root after installing temari_engine, for example:
    python python/temari_engine/benchmarks/bench_batch.py --repeats 5

The script reads only the three published archives from a local directory
(--dist-dir, default dist/ in this checkout) and uses a temporary extraction
cache. No network access or table generation is involved. dist/ is not part of
the public repository: put temari-dataset-v7.0.0.tar.gz,
temari-factors-v2.0.0.tar.gz and temari-factors-ion-v1.0.0.tar.gz there first
(the GitHub Release assets or the Zenodo records; `python -m temari_engine data
info` prints the DOIs).
"""

import argparse
import gc
import json
import os
import platform
import statistics
import sys
import tempfile
import time
import tracemalloc
from pathlib import Path

import numpy as np

from temari_engine import data as td
from temari_engine._factors_reader import _SQUARE_CHUNK, N_NODES


ROOT = Path(__file__).resolve().parents[3]
GIB = 1 << 30


def _ulp_distance(left, right):
    def ordered(value):
        bits = int(np.asarray(value, dtype=np.float64).view(np.uint64))
        return (~bits & ((1 << 64) - 1)) if bits >> 63 else (bits | (1 << 63))
    return abs(ordered(left) - ordered(right))


def _measure(function, repeats, returned_bytes):
    function()  # warm imports, evaluator caches, and allocator pools
    times = []
    for _ in range(repeats):
        gc.collect()
        started = time.perf_counter()
        function()
        times.append(time.perf_counter() - started)
    gc.collect()
    tracemalloc.start()
    result = function()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    del result
    return {"seconds_each": times, "median_seconds": statistics.median(times),
            "tracemalloc_peak_bytes": peak, "returned_bytes": returned_bytes,
            "tracemalloc_scratch_estimate_bytes": max(0, peak - returned_bytes)}


def _scalar_loop(elements, quantity, s_values):
    method_name = quantity

    def run():
        values = np.empty((len(elements), len(s_values)), dtype=np.float64)
        for row, element in enumerate(elements):
            method = getattr(element, method_name)
            for col, s in enumerate(s_values):
                values[row, col] = method(float(s))
        return values
    return run


ARCHIVES = ("temari-dataset-v7.0.0.tar.gz", "temari-factors-v2.0.0.tar.gz", "temari-factors-ion-v1.0.0.tar.gz")


def _benchmark_f(dist_dir, cache_dir, repeats):
    archive = dist_dir / "temari-dataset-v7.0.0.tar.gz"
    dataset = td.open_archive(str(archive), "F", cache_dir=cache_dir)
    channel = dataset.channel(79, "L3")
    e0 = 200.0
    s_values = np.linspace(0.0, channel.s_kin(e0) + 1.0, 1_000_000, dtype=np.float64)
    batch_call = lambda: channel.f_many(e0, s_values)
    probe_indices = np.linspace(0, len(s_values) - 1, 257, dtype=np.int64)
    probe = channel.f_many(e0, s_values[probe_indices])
    max_probe_ulp = max(_ulp_distance(probe.values[i], channel.f(e0, float(s_values[j])).value)
                        for i, j in enumerate(probe_indices))
    if max_probe_ulp > 2:
        raise RuntimeError("F real-data probe differs from scalar by %d ULP (limit 2)" % max_probe_ulp)
    output_bytes = len(s_values) * (8 + 8 + 11 * 4)

    def scalar_call():
        values = np.empty(s_values.shape, dtype=np.float64)
        regions = np.empty(s_values.shape, dtype="U11")
        tail_bounds = np.empty(s_values.shape, dtype=np.float64)
        for i, s in enumerate(s_values):
            result = channel.f(e0, float(s))
            values[i], regions[i], tail_bounds[i] = result.value, result.region, result.tail_bound
        return values, regions, tail_bounds

    batch = _measure(batch_call, repeats, output_bytes)
    scalar = _measure(scalar_call, repeats, output_bytes)
    batch["scratch_upper_bound_bytes"] = 80 * len(s_values)
    batch["regions"] = {str(k): int(v) for k, v in zip(*np.unique(batch_call().regions, return_counts=True))}
    del scalar_call, batch_call
    return {"workload": "F Au L3, E0=200 keV", "points": len(s_values),
            "scalar": scalar, "batch": batch,
            "scalar_sample_max_ulp": max_probe_ulp,
            "speedup": scalar["median_seconds"] / batch["median_seconds"]}


def _benchmark_factors(dist_dir, cache_dir, repeats, series, quantity, ns):
    archive_name = "temari-factors-v2.0.0.tar.gz" if series == "factors" else "temari-factors-ion-v1.0.0.tar.gz"
    dataset = td.open_archive(str(dist_dir / archive_name), series, cache_dir=cache_dir)
    elements = [dataset.element(z, tag) for z, tag in dataset.species()]
    s_values = np.linspace(0.0 if quantity != "fe" else 0.01, 6.0, ns, dtype=np.float64)
    batch_call = lambda: dataset.evaluate_many(elements, quantity, s_values)
    scalar_call = _scalar_loop(elements, quantity, s_values)
    output_bytes = len(elements) * len(s_values) * 8
    result = _measure(batch_call, repeats, output_bytes)
    # 24 per-species vectors of len(s) + the call's cache (6 basis vectors, the squares, and the key and lookup bytes of
    # the query points) + the node arrays and their key and lookup bytes + up to two chunks of Python floats in
    # _square_like_scalar (list slot 8 + float object 24 + array 8 bytes each, rounded up to 48)
    result["scratch_upper_bound_bytes"] = 33 * len(s_values) * 8 + 5 * N_NODES * 8 + 2 * 48 * _SQUARE_CHUNK
    scalar = _measure(scalar_call, repeats, output_bytes)
    del batch_call, scalar_call
    return {"workload": "%s %s, all species" % (series, quantity), "species": len(elements),
            "s_points_per_species": len(s_values), "cells": len(elements) * len(s_values),
            "scalar": scalar, "batch": result,
            "speedup": scalar["median_seconds"] / result["median_seconds"]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--dist-dir", type=Path, default=ROOT / "dist",
                        help="directory holding the three published archives (default: dist/ in this checkout)")
    args = parser.parse_args()
    if args.repeats < 3:
        parser.error("--repeats must be at least 3 to report a meaningful median")
    missing = [name for name in ARCHIVES if not (args.dist_dir / name).is_file()]
    if missing:
        parser.error("archives not found in %s: %s (download the published archives there first; see the module "
                     "docstring)" % (args.dist_dir, ", ".join(missing)))

    with tempfile.TemporaryDirectory(prefix="temari_batch_bench_") as cache_dir:
        results = [
            _benchmark_f(args.dist_dir, cache_dir, args.repeats),
            _benchmark_factors(args.dist_dir, cache_dir, args.repeats, "factors", "fx", 1_000_000 // 86),
            _benchmark_factors(args.dist_dir, cache_dir, args.repeats, "factors-ion", "fe_regular", 1_000_000 // 22),
        ]
    report = {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "cpu_count": os.cpu_count(),
        "repeats": args.repeats,
        "acceptance": {"minimum_speedup": 10.0, "maximum_scratch_bytes": GIB},
        "results": results,
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    failed = any(item["speedup"] < 10.0 or item["batch"]["scratch_upper_bound_bytes"] > GIB for item in results)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
