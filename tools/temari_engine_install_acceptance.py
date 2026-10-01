#!/usr/bin/env python3
"""Install acceptance for the published temari_engine wheel.

Run this from an isolated Python environment outside a Temari checkout.  The
first invocation fills a new cache.  A second process runs with
TEMARI_OFFLINE=1 and compares all stable results with the first invocation.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import struct
import sys
import time
import traceback

import temari_engine as te
from temari_engine import data as td


BATCH_POINTS = 2000

EXPECTED = {
    "temari_engine": "0.5.0",
    "F": {
        "version": "7.0.0",
        "doi": "10.5281/zenodo.22643468",
        "files": 525,
        "Au_L3_200keV_s1": 0.862270670072035,
    },
    "factors": {
        "version": "2.0.0",
        "doi": "10.5281/zenodo.22820415",
        "files": 86,
        "Fe_fx_0.5": 11.510029462,
        "Fe_fe_0.5": 1.3871921689,
    },
    "factors-ion": {
        "version": "1.0.0",
        "doi": "10.5281/zenodo.22820492",
        "files": 22,
        "O2-_first": ["29d4d8d15a61afb0", 4, 2.2992044845],
    },
    "Au_LA": {
        "sigma_nm2": 1.0474175553264528e-08,
        # SHA-256 of all 321 values formatted as ".12e", one per LF line.
        "shape_sha256_13_digits": "5148dd642814b7c5214229238c2c95977583d3f58463a443669400ffbf1df73a",
        "weights": {
            "L1": 0.15501451976179687,
            "L2": 0.040825352441135575,
            "L3": 0.8041601277970676,
        },
        "shape_samples": {
            0: 1.0,
            1: 0.9996109742755512,
            80: 0.22347918697008423,
            160: 0.009846315655017326,
            240: 0.002416167166769212,
            320: 0.0009600236020084481,
        },
    },
}


def _close(actual, expected):
    return math.isfinite(actual) and math.isclose(
        actual, expected, rel_tol=1e-12, abs_tol=1e-15
    )


def _shape_digest(values):
    normalized = ("\n".join(format(value, ".12e") for value in values) + "\n").encode("ascii")
    return hashlib.sha256(normalized).hexdigest()


def _is_below(path, parent):
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def _ulp(left, right):
    def ordered(value):
        bits = struct.unpack(">Q", struct.pack(">d", float(value)))[0]
        return (~bits & ((1 << 64) - 1)) if bits >> 63 else (bits | (1 << 63))
    return abs(ordered(left) - ordered(right))


def _batch_measure(au_l3, factor_sets):
    """0.5.0: the NumPy batch APIs on the published tables against the scalar methods on this platform.

    F batch values are promised within 2 ULP of scalar (same regions); dataset-factors batch values are bit-identical."""
    import numpy as np

    s_values = np.random.default_rng(20261001).uniform(0.0, 6.0, BATCH_POINTS)
    many = au_l3.f_many(200.0, s_values)
    scalar = [au_l3.f(200.0, s) for s in s_values.tolist()]
    cells = mismatches = 0
    species = []
    for dataset in factor_sets:
        elements = [dataset.element(z, tag) for z, tag in dataset.species()]
        species.append(len(elements))
        for quantity in ("fx", "fe", "fe_regular"):
            rows = dataset.evaluate_many(elements, quantity, s_values).values.tolist()
            for row, element in zip(rows, elements):
                method = getattr(element, quantity)
                for actual, s in zip(row, s_values.tolist()):
                    cells += 1
                    mismatches += struct.pack("<d", actual) != struct.pack("<d", method(s))
    return {
        "points": len(scalar),
        "F_Au_L3_200keV_max_ulp": max(_ulp(a, b.value) for a, b in zip(many.values.tolist(), scalar)),
        "F_Au_L3_200keV_regions_equal": list(many.regions) == [b.region for b in scalar],
        "factors_species": species,
        "factors_cells": cells,
        "factors_bit_mismatches": mismatches,
    }


def _measure():
    f = td.open("F", verify="full")
    au_l3 = f.channel(79, "L3")
    r = au_l3.f(200.0, 1.0)
    f_result = {
        "files_verified": len(f._store.entries),
        "doi": f.source.doi,
        "version": f.source.version,
        "verified": f.source.verified,
        "Au_L3_200keV_s1": [r.value, r.region, r.tail_bound],
        "Au_L3_200keV_s20": [au_l3.f(200.0, 20.0).region],
        "citation": f.source.citation(),
    }

    factors = td.open("factors", verify="full")
    fe = factors.element(26)
    factors_result = {
        "files_verified": len(factors._store.entries),
        "doi": factors.source.doi,
        "version": factors.source.version,
        "verified": factors.source.verified,
        "Fe_fx_0.5": fe.fx(0.5),
        "Fe_fe_0.5": fe.fe(0.5),
    }

    ions = td.open("factors-ion", verify="full")
    candidates = ions.find(8, charge=-2)
    try:
        ions.find_one(8, -2)
        ambiguous = False
    except te.SpeciesLookupError as exc:
        ambiguous = len(exc.candidates)
    first = candidates[0]
    ions_result = {
        "files_verified": len(ions._store.entries),
        "doi": ions.source.doi,
        "version": ions.source.version,
        "verified": ions.source.verified,
        "O2-_candidates": len(candidates),
        "find_one_refused_with": ambiguous,
        "O2-_first": [
            first.species_tag,
            first.coordination_number,
            ions.element(8, first.species_tag).fx(0.5),
        ],
    }

    from temari_engine import xray_lines

    line_response = xray_lines.line_response(79, 200.0, f.path)
    la = line_response.series("LA")
    line_result = {
        "sigma_nm2": la.sigma_nm2,
        "weights": la.weights,
        "shape_points": len(la.shape),
        "shape": la.shape,
        "shape_sha256_13_digits": _shape_digest(la.shape),
    }
    return {
        "F": f_result,
        "factors": factors_result,
        "factors-ion": ions_result,
        "Au_LA": line_result,
        "batch": _batch_measure(au_l3, (factors, ions)),
    }


def _data_checks(data):
    f = data["F"]
    factors = data["factors"]
    ions = data["factors-ion"]
    la = data["Au_LA"]
    batch = data["batch"]
    return {
        "batch_api": (
            batch["points"] == BATCH_POINTS
            and batch["F_Au_L3_200keV_max_ulp"] <= 2
            and batch["F_Au_L3_200keV_regions_equal"] is True
            and batch["factors_species"] == [86, 22]
            and batch["factors_cells"] == BATCH_POINTS * 3 * (86 + 22)
            and batch["factors_bit_mismatches"] == 0
        ),
        "F_identity": (
            f["files_verified"] == EXPECTED["F"]["files"]
            and f["doi"] == EXPECTED["F"]["doi"]
            and f["version"] == EXPECTED["F"]["version"]
            and f["verified"] == "archive+manifest"
        ),
        "F_representative_values": (
            _close(f["Au_L3_200keV_s1"][0], EXPECTED["F"]["Au_L3_200keV_s1"])
            and f["Au_L3_200keV_s1"][1:] == ["tabulated", 0.0]
            and f["Au_L3_200keV_s20"] == ["unrecorded"]
        ),
        "factors_identity": (
            factors["files_verified"] == EXPECTED["factors"]["files"]
            and factors["doi"] == EXPECTED["factors"]["doi"]
            and factors["version"] == EXPECTED["factors"]["version"]
            and factors["verified"] == "archive+manifest"
        ),
        "factors_representative_values": (
            _close(factors["Fe_fx_0.5"], EXPECTED["factors"]["Fe_fx_0.5"])
            and _close(factors["Fe_fe_0.5"], EXPECTED["factors"]["Fe_fe_0.5"])
        ),
        "ion_identity": (
            ions["files_verified"] == EXPECTED["factors-ion"]["files"]
            and ions["doi"] == EXPECTED["factors-ion"]["doi"]
            and ions["version"] == EXPECTED["factors-ion"]["version"]
            and ions["verified"] == "archive+manifest"
        ),
        "ion_lookup": (
            ions["O2-_candidates"] == 5
            and ions["find_one_refused_with"] == 5
            and ions["O2-_first"][:2] == EXPECTED["factors-ion"]["O2-_first"][:2]
            and _close(ions["O2-_first"][2], EXPECTED["factors-ion"]["O2-_first"][2])
        ),
        "Au_LA": (
            _close(la["sigma_nm2"], EXPECTED["Au_LA"]["sigma_nm2"])
            and la["shape_points"] == 321
            and set(la["weights"]) == {"L1", "L2", "L3"}
            and all(
                _close(la["weights"][shell], expected)
                for shell, expected in EXPECTED["Au_LA"]["weights"].items()
            )
            and len(la["shape"]) == 321
            and all(math.isfinite(value) for value in la["shape"])
            and _shape_digest(la["shape"]) == EXPECTED["Au_LA"]["shape_sha256_13_digits"]
            and all(
                _close(la["shape"][index], expected)
                for index, expected in EXPECTED["Au_LA"]["shape_samples"].items()
            )
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--forbid-path", required=True)
    parser.add_argument("--expect-offline", action="store_true")
    parser.add_argument("--compare-with")
    args = parser.parse_args()

    started = time.time()
    module_file = Path(te.__file__).resolve()
    forbidden = Path(args.forbid_path).resolve()
    cache = Path(td.cache_root()).resolve()
    expected_cache = os.environ.get("TEMARI_CACHE_DIR")
    offline_env = os.environ.get("TEMARI_OFFLINE", "")
    result = {
        "environment": {
            "python": sys.version.split()[0],
            "platform": sys.platform,
            "temari_engine": te.__version__,
            "module_file": str(module_file),
            "forbidden_path": str(forbidden),
            "cache": str(cache),
            "offline_env": offline_env,
            "julia": shutil.which("julia"),
        },
        "data": None,
        "checks": {
            "package_version": te.__version__ == EXPECTED["temari_engine"],
            "linux_platform": sys.platform == "linux",
            "module_outside_repo": not _is_below(module_file, forbidden),
            "installed_in_site_packages": "site-packages" in module_file.parts,
            "no_development_pythonpath": not os.environ.get("PYTHONPATH"),
            "julia_absent": shutil.which("julia") is None,
            "offline_mode": offline_env == ("1" if args.expect_offline else ""),
            "dedicated_cache": bool(expected_cache) and cache == Path(expected_cache).resolve(),
        },
        "error": None,
    }

    try:
        result["data"] = _measure()
        result["checks"].update(_data_checks(result["data"]))
        if args.compare_with:
            with open(args.compare_with, "r", encoding="utf-8") as stream:
                previous = json.load(stream)
            result["checks"]["online_offline_equal"] = result["data"] == previous.get("data")
    except Exception as exc:  # Keep a machine-readable record for a failed Actions run.
        result["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }

    result["elapsed_s"] = round(time.time() - started, 1)
    result["accept"] = result["error"] is None and all(result["checks"].values())
    rendered = json.dumps(result, indent=1, sort_keys=True, ensure_ascii=False) + "\n"
    Path(args.output).write_text(rendered, encoding="utf-8", newline="\n")
    print(rendered, end="")
    print("ACCEPT=%s" % ("pass" if result["accept"] else "fail"))
    return 0 if result["accept"] else 1


if __name__ == "__main__":
    sys.exit(main())
