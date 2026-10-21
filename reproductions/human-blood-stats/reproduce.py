#!/usr/bin/env python3
"""Count repository sample/factor rows in studies tagged as human.

Counts, per repository:

* the number of studies tagged with a human organism/species,
* the total number of returned sample/factor rows in those studies,
* the number of rows whose source/tissue text matches the broad blood rule.

The study-level tag does not establish that every row is a human biological
specimen. Counts are rows, not unique people.

Method
------
* MetaboLights: ``search(organisms=["Homo sapiens"])`` then a deep ``inspect()``
  of every human-tagged study; every parsed ISA sample-sheet row is counted,
  and ``organism_part`` is used for the broad blood-label tally. Studies that fail to
  inspect are recorded with their search-index ``sample_count`` so the
  undercount is explicit.
* Metabolomics Workbench: the whole-index ``/ST/summary`` corpus is filtered to
  ``species == "Homo sapiens"``; every selected study's ``factors`` endpoint is
  fetched, every returned row is counted, and ``sample_source`` is used for the
  broad blood-label tally.

Blood classification
--------------------
A sample is "blood of some sort" when its tissue/source contains
``blood``/``plasma``/``serum``. ``seminal plasma`` is explicitly excluded (it
matches on "plasma" but is not blood). "canonical" means exactly
``blood plasma`` / ``blood serum``.

Usage
-----
    # after installing the root requirements.txt:
    python reproductions/human-blood-stats/reproduce.py
    python reproductions/human-blood-stats/reproduce.py --workers 8 --out out

The script is read-only: it never writes into the metabo-search repository and
downloads nothing except the public metadata it counts.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

from metabo_search import MetaboLights, Workbench

#: Commit of the metabo-search library these numbers were produced with.
METABO_SEARCH_COMMIT = "4ee8d01525a898b83dfcd3ec408e372f013236f8"

#: Substrings that would otherwise match "plasma" but are not blood.
_BLOOD_EXCLUDE = ("seminal",)

_CANONICAL = frozenset({"blood plasma", "blood serum"})


def is_blood(term: str) -> bool:
    """True when tissue/source text matches the broad blood-label rule."""
    text = (term or "").strip().lower()
    if not text or any(bad in text for bad in _BLOOD_EXCLUDE):
        return False
    return "blood" in text or "plasma" in text or "serum" in text or text in ("plasma", "serum")


def is_canonical(term: str) -> bool:
    """True for the canonical MetaboLights biofluids ``blood plasma``/``blood serum``."""
    return (term or "").strip().lower() in _CANONICAL


# --------------------------------------------------------------------------
# MetaboLights
# --------------------------------------------------------------------------


def _inspect_one(ml: MetaboLights, hit: Any) -> dict[str, Any]:
    """Inspect a single study, never raising: failures are recorded, not lost."""
    base = {
        "study_id": hit.study_id,
        "index_samples": hit.sample_count or 0,
        "parsed_samples": 0,
        "error": None,
        "study_parts": [t.term for t in hit.organism_parts],
        "parts": {},
    }
    try:
        detail = ml.inspect([hit], workers=1)[0]
    except Exception as exc:  # noqa: BLE001 - report per study, keep the batch
        base["error"] = f"{type(exc).__name__}: {exc}"
        return base
    if detail.error:
        # Best-effort inspect flags studies with no readable metadata instead
        # of raising (e.g. a study whose mirror has no ISA files).
        base["error"] = detail.error
        return base
    base["parsed_samples"] = len(detail.samples)
    base["parts"] = dict(Counter(s.organism_part for s in detail.samples))
    return base


def collect_metabolights(workers: int, limit: int) -> dict[str, Any]:
    """Count rows in human-tagged MetaboLights studies."""
    ml = MetaboLights()
    try:
        hits = ml.search("", organisms=["Homo sapiens"], limit=limit)
        with ThreadPoolExecutor(max_workers=workers) as pool:
            per_study = list(pool.map(lambda h: _inspect_one(ml, h), hits))
    finally:
        ml.close()

    # A study is "unreliable" when inspect raised, or when it silently parsed
    # zero samples although the search index reports some (the empty-metadata
    # case: ISA files missing from the mirror). Both are undercounts.
    reliable = [
        r for r in per_study if r["error"] is None and not (r["index_samples"] > 0 and r["parsed_samples"] == 0)
    ]
    unreliable = [r for r in per_study if r not in reliable]

    total_samples = sum(r["parsed_samples"] for r in reliable)
    total_estimate = total_samples + sum(r["index_samples"] for r in unreliable)
    blood = canonical = 0
    for record in reliable:
        for term, count in record["parts"].items():
            if is_blood(term):
                blood += count
            if is_canonical(term):
                canonical += count
    # For unreliable studies, fall back to the study-level facet from the
    # search index (cannot see individual samples).
    blood_estimate = blood
    canonical_estimate = canonical
    for record in unreliable:
        parts = record["study_parts"]
        if any(is_blood(t) for t in parts):
            blood_estimate += record["index_samples"]
        if any(is_canonical(t) for t in parts):
            canonical_estimate += record["index_samples"]

    return {
        "repository": "metabolights",
        "human_studies": len(hits),
        "human_samples_total": total_samples,
        "human_samples_total_estimate": total_estimate,
        "human_samples_blood": blood,
        "human_samples_blood_estimate": blood_estimate,
        "human_samples_canonical": canonical,
        "human_samples_canonical_estimate": canonical_estimate,
        "studies_inspected_ok": len(reliable),
        "studies_unreliable": len(unreliable),
        "unreliable_studies": [
            {
                "study_id": r["study_id"],
                "index_samples": r["index_samples"],
                "parsed_samples": r["parsed_samples"],
                "study_parts": r["study_parts"],
                "error": r["error"],
            }
            for r in unreliable
        ],
        "per_study": per_study,
    }


# --------------------------------------------------------------------------
# Metabolomics Workbench
# --------------------------------------------------------------------------


def _factors_one(wb: Workbench, summary: Any) -> dict[str, Any]:
    """Fetch one study's factor rows, never raising."""
    base = {
        "study_id": summary.study_id,
        "summary_samples": summary.number_of_samples or 0,
        "factor_rows": 0,
        "error": None,
        "sources": {},
    }
    try:
        rows = wb.samples(summary.study_id)
    except Exception as exc:  # noqa: BLE001 - report per study, keep the batch
        base["error"] = f"{type(exc).__name__}: {exc}"
        return base
    base["factor_rows"] = len(rows)
    base["sources"] = dict(Counter(r.sample_source for r in rows))
    return base


def collect_workbench(workers: int) -> dict[str, Any]:
    """Count rows in human-tagged Workbench studies."""
    wb = Workbench()
    try:
        summaries = wb.all_studies()
        human = [s for s in summaries.values() if s.species == "Homo sapiens"]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            per_study = list(pool.map(lambda s: _factors_one(wb, s), human))
    finally:
        wb.close()

    total_samples = sum(r["factor_rows"] for r in per_study)
    blood = 0
    for record in per_study:
        for source, count in record["sources"].items():
            if is_blood(source):
                blood += count
    # Studies whose factors endpoint is empty: cover them with the summary count.
    missing = sum(r["summary_samples"] for r in per_study if r["factor_rows"] == 0)

    return {
        "repository": "metabolomics_workbench",
        "human_studies": len(human),
        "human_samples_total": total_samples,
        "human_samples_total_estimate": total_samples + missing,
        "human_samples_blood": blood,
        "human_samples_blood_estimate": blood,
        "studies_failed": sum(1 for r in per_study if r["error"] is not None),
        "per_study": per_study,
    }


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------


def _fmt(value: int) -> str:
    return f"{value:,}"


def _print_table(ml: dict[str, Any], wb: dict[str, Any]) -> None:
    header = (
        f"{'repository':<26}{'human-tagged':>14}{'parsed rows':>16}"
        f"{'blood-label rows':>18}{'est. rows':>14}{'est. blood':>12}"
    )
    print(header)
    print("-" * len(header))
    for block, label in ((ml, "MetaboLights"), (wb, "Metabolomics Workbench")):
        print(
            f"{label:<26}{_fmt(block['human_studies']):>14}"
            f"{_fmt(block['human_samples_total']):>16}"
            f"{_fmt(block['human_samples_blood']):>18}"
            f"{_fmt(block['human_samples_total_estimate']):>14}"
            f"{_fmt(block['human_samples_blood_estimate']):>12}"
        )
    total_studies = ml["human_studies"] + wb["human_studies"]
    total_samples = ml["human_samples_total"] + wb["human_samples_total"]
    total_blood = ml["human_samples_blood"] + wb["human_samples_blood"]
    est_samples = ml["human_samples_total_estimate"] + wb["human_samples_total_estimate"]
    est_blood = ml["human_samples_blood_estimate"] + wb["human_samples_blood_estimate"]
    print("-" * len(header))
    print(
        f"{'COMBINED':<26}{_fmt(total_studies):>14}{_fmt(total_samples):>16}"
        f"{_fmt(total_blood):>18}{_fmt(est_samples):>14}{_fmt(est_blood):>12}"
    )
    print(
        "\n'est.' adds study-index row counts for studies whose sample metadata\n"
        "could not be read (see studies_unreliable). These are approximate row\n"
        "totals, not sample-by-sample classifications or unique participants."
    )
    if est_samples:
        print(f"\nblood-label share of estimated rows: {est_blood / est_samples:.1%}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="parallel study fetches (4 is the measured sweet spot)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=100_000,
        help="max MetaboLights search hits to fetch (default: all)",
    )
    parser.add_argument("--out", type=Path, default=Path("out"), help="directory for JSON results")
    parser.add_argument("--skip-metabolights", action="store_true")
    parser.add_argument("--skip-workbench", action="store_true")
    args = parser.parse_args()

    if args.skip_metabolights and args.skip_workbench:
        parser.error("nothing to do: both repositories skipped")

    print(f"metabo-search commit: {METABO_SEARCH_COMMIT}\n")

    ml: dict[str, Any] = {
        "repository": "metabolights",
        "human_studies": 0,
        "human_samples_total": 0,
        "human_samples_total_estimate": 0,
        "human_samples_blood": 0,
        "human_samples_blood_estimate": 0,
    }
    wb: dict[str, Any] = {
        "repository": "metabolomics_workbench",
        "human_studies": 0,
        "human_samples_total": 0,
        "human_samples_total_estimate": 0,
        "human_samples_blood": 0,
        "human_samples_blood_estimate": 0,
    }

    args.out.mkdir(parents=True, exist_ok=True)

    if not args.skip_metabolights:
        print("== MetaboLights ==", flush=True)
        ml = collect_metabolights(args.workers, args.limit)
        (args.out / "metabolights.json").write_text(json.dumps(ml, indent=2))
        print(
            f"  {_fmt(ml['human_studies'])} human-tagged studies; "
            f"{_fmt(ml['human_samples_total'])} rows; "
            f"{_fmt(ml['human_samples_blood'])} blood-label rows "
            f"({ml['studies_unreliable']} studies with unreadable metadata; "
            f"est. {_fmt(ml['human_samples_total_estimate'])} rows, "
            f"{_fmt(ml['human_samples_blood_estimate'])} blood-label rows)\n",
            flush=True,
        )

    if not args.skip_workbench:
        print("== Metabolomics Workbench ==", flush=True)
        wb = collect_workbench(args.workers)
        (args.out / "workbench.json").write_text(json.dumps(wb, indent=2))
        print(
            f"  {_fmt(wb['human_studies'])} human-tagged studies; "
            f"{_fmt(wb['human_samples_total'])} rows; "
            f"{_fmt(wb['human_samples_blood'])} blood-label rows\n",
            flush=True,
        )

    _print_table(ml, wb)
    summary = {
        "metabo_search_commit": METABO_SEARCH_COMMIT,
        "metabolights": {k: v for k, v in ml.items() if k != "per_study"},
        "workbench": {k: v for k, v in wb.items() if k != "per_study"},
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nJSON written to {args.out}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
