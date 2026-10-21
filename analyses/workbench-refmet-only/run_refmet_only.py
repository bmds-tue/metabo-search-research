#!/usr/bin/env python3
"""End-to-end, Workbench-only human-study RefMet-name overlap analysis.

Default mode starts from the live Workbench API, verifies human samples, fetches
only RefMet names, and writes a DuckDB snapshot, CSVs, manifest, and SVG plot.
Use ``--recompute RUN_ID`` to rebuild scores/exports from an existing snapshot
without network access.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import re
import sys
import time
import unicodedata
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import duckdb
import pyarrow as pa

ROOT = Path(__file__).resolve().parent
VENDORED_SRC = ROOT / "vendor" / "metabo-search" / "src"
sys.path.insert(0, str(VENDORED_SRC))
from metabo_search import Workbench  # noqa: E402
from metabo_search.workbench.models import StudyDetail, StudySummary  # noqa: E402

DB_DEFAULT = ROOT / "workbench_refmet_only.duckdb"
MIN_REFMET_NAMES = 100
WORKERS = 4
BATCH_SIZE = 40
RUN_ID_FORMAT = "%Y%m%dT%H%M%SZ"

SCHEMA = """
CREATE TABLE IF NOT EXISTS run_info (
  run_id VARCHAR PRIMARY KEY,
  started_at VARCHAR,
  completed_at VARCHAR,
  status VARCHAR,
  metabo_search_version VARCHAR,
  workbench_summary_count INTEGER,
  human_candidate_count INTEGER,
  human_sample_study_count INTEGER,
  eligible_study_count INTEGER,
  pair_count BIGINT,
  mean_jaccard DOUBLE,
  population_variance DOUBLE,
  sample_variance DOUBLE,
  errors_json VARCHAR,
  manifest_json VARCHAR
);
CREATE TABLE IF NOT EXISTS study (
  run_id VARCHAR,
  study_id VARCHAR,
  title VARCHAR,
  species VARCHAR,
  sample_status VARCHAR,
  sample_count INTEGER,
  human_sample_count INTEGER,
  human_detection_method VARCHAR,
  sample_error VARCHAR,
  metabolite_status VARCHAR,
  refmet_row_count INTEGER,
  unique_refmet_name_count INTEGER,
  qualifies BOOLEAN,
  metabolite_error VARCHAR,
  PRIMARY KEY (run_id, study_id)
);
CREATE TABLE IF NOT EXISTS human_sample (
  run_id VARCHAR,
  study_id VARCHAR,
  sample_ordinal INTEGER,
  sample_name VARCHAR,
  organism_evidence VARCHAR,
  sample_source VARCHAR,
  evidence VARCHAR,
  detection_method VARCHAR,
  PRIMARY KEY (run_id, study_id, sample_ordinal)
);
CREATE TABLE IF NOT EXISTS study_refmet_name (
  run_id VARCHAR,
  study_id VARCHAR,
  name_key VARCHAR,
  refmet_name VARCHAR,
  PRIMARY KEY (run_id, study_id, name_key)
);
CREATE TABLE IF NOT EXISTS pairwise_jaccard (
  run_id VARCHAR,
  study1_id VARCHAR,
  study2_id VARCHAR,
  intersection_count INTEGER,
  union_count INTEGER,
  jaccard DOUBLE,
  PRIMARY KEY (run_id, study1_id, study2_id)
);
CREATE TABLE IF NOT EXISTS summary (
  run_id VARCHAR PRIMARY KEY,
  human_sample_study_count INTEGER,
  studies_with_refmet_names INTEGER,
  eligible_study_count INTEGER,
  pair_count BIGINT,
  mean_jaccard DOUBLE,
  population_variance DOUBLE,
  sample_variance DOUBLE,
  zero_overlap_pairs BIGINT,
  median_jaccard DOUBLE,
  minimum_jaccard DOUBLE,
  maximum_jaccard DOUBLE
);
"""

PAIR_INSERT = """
WITH eligible AS (
  SELECT study_id FROM study WHERE run_id=? AND qualifies
),
names AS (
  SELECT n.study_id,n.name_key FROM study_refmet_name n
  JOIN eligible e USING (study_id) WHERE n.run_id=?
),
sizes AS (
  SELECT study_id,count(*)::INTEGER AS n FROM names GROUP BY study_id
),
intersections AS (
  SELECT a.study_id AS id1,b.study_id AS id2,count(*)::INTEGER AS n
  FROM names a JOIN names b ON a.name_key=b.name_key AND a.study_id<b.study_id
  GROUP BY a.study_id,b.study_id
)
INSERT INTO pairwise_jaccard
  (run_id,study1_id,study2_id,intersection_count,union_count,jaccard)
SELECT ?,a.study_id,b.study_id,coalesce(i.n,0),sa.n+sb.n-coalesce(i.n,0),
       coalesce(i.n,0)::DOUBLE/(sa.n+sb.n-coalesce(i.n,0))
FROM eligible a CROSS JOIN eligible b
JOIN sizes sa ON sa.study_id=a.study_id
JOIN sizes sb ON sb.study_id=b.study_id
LEFT JOIN intersections i ON i.id1=a.study_id AND i.id2=b.study_id
WHERE a.study_id<b.study_id
ON CONFLICT (run_id,study1_id,study2_id) DO UPDATE SET
 intersection_count=excluded.intersection_count,
 union_count=excluded.union_count,
 jaccard=excluded.jaccard
"""

_HUMAN = re.compile(r"\bhomo\s+sapiens\b|\bh\.?\s*sapiens\b|\bhuman\b", re.IGNORECASE)
_TAXON_LABEL = re.compile(r"species|organism|taxon", re.IGNORECASE)
_UNKNOWN = {"", "-", "unknown", "not applicable", "n/a"}


def human_label(value: str | None) -> bool:
    return bool(_HUMAN.search(value or ""))


def single_human_summary(species: str) -> bool:
    if not human_label(species):
        return False
    labels = [part.strip() for part in re.split(r"[,;/|]", species) if part.strip()]
    return bool(labels) and all(human_label(part) for part in labels)


def canonical_name(value: str) -> tuple[str, str] | None:
    display = " ".join(unicodedata.normalize("NFKC", value or "").split())
    return (display.casefold(), display) if display else None


def find_human_samples(
    detail: StudyDetail, summary: StudySummary
) -> tuple[list[tuple[int, str, str, str, str, str]], str]:
    """Classify human samples; explicit non-human taxonomy overrides summary fallback."""
    found: list[tuple[int, str, str, str, str, str]] = []
    fallback = single_human_summary(summary.species)
    for ordinal, sample in enumerate(detail.samples):
        taxon_values = [
            (label, value.strip())
            for label, value in sample.factor_values.items()
            if _TAXON_LABEL.search(label)
            and value.strip().casefold() not in _UNKNOWN
        ]
        human_taxa = [(label, value) for label, value in taxon_values if human_label(value)]
        if human_taxa:
            evidence = "; ".join(f"{label}={value}" for label, value in human_taxa)
            found.append((ordinal, sample.local_sample_id or sample.mb_sample_id,
                          evidence, sample.sample_source, evidence, "sample_level"))
        elif taxon_values:
            continue
        elif human_label(sample.sample_source):
            evidence = f"sample_source={sample.sample_source}"
            found.append((ordinal, sample.local_sample_id or sample.mb_sample_id,
                          sample.sample_source, sample.sample_source, evidence, "sample_level"))
        elif fallback:
            evidence = f"summary species={summary.species}; sample taxonomy absent"
            found.append((ordinal, sample.local_sample_id or sample.mb_sample_id,
                          summary.species, sample.sample_source, evidence,
                          "single_species_summary_fallback"))
    methods = sorted({row[5] for row in found})
    return found, "+".join(methods) if methods else "none"


def to_arrow_insert(
    con: duckdb.DuckDBPyConnection,
    table: str,
    columns: tuple[str, ...],
    rows: list[tuple[Any, ...]],
) -> None:
    if not rows:
        return
    arrow = pa.Table.from_pylist([dict(zip(columns, row)) for row in rows])
    con.register("_incoming", arrow)
    try:
        con.execute(f"INSERT INTO {table} ({','.join(columns)}) SELECT * FROM _incoming")
    finally:
        con.unregister("_incoming")


def retry_failed_samples(wb: Workbench, ids: list[str]) -> list[StudyDetail]:
    if not ids:
        return []
    time.sleep(1)
    return list(wb.inspect(ids, workers=1, include_summary=False, include_analyses=False,
                           include_metabolites=False, on_error="collect"))


def retry_failed_metabolites(wb: Workbench, ids: list[str]) -> list[StudyDetail]:
    if not ids:
        return []
    time.sleep(1)
    return list(wb.inspect(ids, workers=1, include_summary=False, include_samples=False,
                           include_analyses=False, include_metabolites=True,
                           on_error="collect"))


def make_pairs(con: duckdb.DuckDBPyConnection, run_id: str) -> dict[str, Any]:
    con.execute("DELETE FROM pairwise_jaccard WHERE run_id=?", [run_id])
    con.execute(PAIR_INSERT, [run_id, run_id, run_id])
    metrics = con.execute(
        "SELECT count(*),avg(jaccard),var_pop(jaccard),var_samp(jaccard), "
        "count(*) FILTER (WHERE jaccard=0),median(jaccard),min(jaccard),max(jaccard) "
        "FROM pairwise_jaccard WHERE run_id=?", [run_id]
    ).fetchone()
    pairs, mean, var_pop, var_samp, zero_pairs, median, minimum, maximum = metrics
    eligible = con.execute(
        "SELECT count(*) FROM study WHERE run_id=? AND qualifies", [run_id]
    ).fetchone()[0]
    confirmed = con.execute(
        "SELECT count(*) FROM study WHERE run_id=? AND sample_status='success' AND human_sample_count>0",
        [run_id],
    ).fetchone()[0]
    refmet_studies = con.execute(
        "SELECT count(*) FROM study WHERE run_id=? AND unique_refmet_name_count>0", [run_id]
    ).fetchone()[0]
    con.execute("""
      INSERT INTO summary VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
      ON CONFLICT (run_id) DO UPDATE SET
        human_sample_study_count=excluded.human_sample_study_count,
        studies_with_refmet_names=excluded.studies_with_refmet_names,
        eligible_study_count=excluded.eligible_study_count,pair_count=excluded.pair_count,
        mean_jaccard=excluded.mean_jaccard,population_variance=excluded.population_variance,
        sample_variance=excluded.sample_variance,zero_overlap_pairs=excluded.zero_overlap_pairs,
        median_jaccard=excluded.median_jaccard,minimum_jaccard=excluded.minimum_jaccard,
        maximum_jaccard=excluded.maximum_jaccard
    """, [run_id, confirmed, refmet_studies, eligible, pairs, mean, var_pop, var_samp,
          zero_pairs, median, minimum, maximum])
    return {
        "human_sample_study_count": confirmed,
        "studies_with_refmet_names": refmet_studies,
        "eligible_study_count": eligible,
        "pair_count": pairs,
        "mean_jaccard": mean,
        "population_variance": var_pop,
        "sample_variance": var_samp,
        "zero_overlap_pairs": zero_pairs,
        "median_jaccard": median,
        "minimum_jaccard": minimum,
        "maximum_jaccard": maximum,
    }


def write_csv(
    con: duckdb.DuckDBPyConnection, query: str, params: list[Any], path: Path,
    header: list[str],
) -> None:
    rows = con.execute(query, params).fetchall()
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def make_plot(con: duckdb.DuckDBPyConnection, run_id: str, stats: dict[str, Any]) -> None:
    scores = [r[0] for r in con.execute(
        "SELECT jaccard FROM pairwise_jaccard WHERE run_id=? ORDER BY jaccard", [run_id]
    ).fetchall()]
    n = len(scores)
    if not n:
        return
    width, height = 1100, 560
    top, plot_h = 104, 340
    panel_w, gap, left1 = 455, 70, 85
    left2 = left1 + panel_w + gap
    bins_full, bins_zoom = 50, 25
    zoom_max = 0.25
    counts_full = [0] * bins_full
    counts_zoom = [0] * bins_zoom
    for value in scores:
        index = min(int(value * bins_full), bins_full - 1)
        counts_full[index] += 1
        if value <= zoom_max:
            counts_zoom[min(int(value / zoom_max * bins_zoom), bins_zoom - 1)] += 1
    ymax = max(counts_full) / n * 100
    ymax_zoom = max(counts_zoom) / n * 100
    mean, median = stats["mean_jaccard"], stats["median_jaccard"]
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#243447}.title{font-size:22px;font-weight:bold}.subtitle{font-size:14px}.axis{stroke:#34495e;stroke-width:1}.grid{stroke:#dfe6e9;stroke-width:1}.bar{fill:#4386b5}.mean{stroke:#c0392b;stroke-width:2}.median{stroke:#16845b;stroke-width:2;stroke-dasharray:5 4}.tick{font-size:12px}.label{font-size:14px}</style>',
             f'<text x="{width/2}" y="32" text-anchor="middle" class="title">Workbench human-study Jaccard distribution (RefMet names only)</text>',
             f'<text x="{width/2}" y="57" text-anchor="middle" class="subtitle">{stats["eligible_study_count"]:,} studies · {n:,} unordered pairs · mean {mean:.4f} · median {median:.4f}</text>']

    def panel(left: int, x_max: float, counts: list[int], bins: int, max_percent: float, subtitle: str, zoom: bool) -> None:
        bottom = top + plot_h
        parts.append(f'<text x="{left+panel_w/2}" y="88" text-anchor="middle" class="label">{subtitle}</text>')
        for tick in range(5):
            y = bottom - plot_h * tick / 4
            pct = max_percent * tick / 4
            parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left+panel_w}" y2="{y:.1f}" class="grid"/>')
            parts.append(f'<text x="{left-10}" y="{y+4:.1f}" text-anchor="end" class="tick">{pct:.1f}%</text>')
        bw = panel_w / bins
        for i, count in enumerate(counts):
            bar_h = (count / n * 100 / max_percent * plot_h) if max_percent else 0
            x = left + i * bw
            y = bottom - bar_h
            parts.append(f'<rect x="{x+0.4:.2f}" y="{y:.2f}" width="{max(0,bw-0.8):.2f}" height="{bar_h:.2f}" class="bar"/>')
        parts.append(f'<line x1="{left}" y1="{bottom}" x2="{left+panel_w}" y2="{bottom}" class="axis"/>')
        parts.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" class="axis"/>')
        x_ticks = [0,.2,.4,.6,.8,1] if not zoom else [0,.05,.10,.15,.20,.25]
        for tick in x_ticks:
            x = left + tick / x_max * panel_w
            parts.append(f'<line x1="{x:.1f}" y1="{bottom}" x2="{x:.1f}" y2="{bottom+5}" class="axis"/>')
            parts.append(f'<text x="{x:.1f}" y="{bottom+22}" text-anchor="middle" class="tick">{tick:.2f}</text>')
        parts.append(f'<text x="{left+panel_w/2}" y="{bottom+47}" text-anchor="middle" class="label">Jaccard index</text>')
        for value, klass in ((mean,"mean"),(median,"median")):
            if value <= x_max:
                x = left + value / x_max * panel_w
                parts.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{bottom}" class="{klass}"/>')

    panel(left1,1.0,counts_full,bins_full,ymax*1.08,"Full range (0–1)",False)
    panel(left2,zoom_max,counts_zoom,bins_zoom,ymax_zoom*1.08,"Zoom near zero (0–0.25)",True)
    parts.extend([
        f'<line x1="355" y1="520" x2="382" y2="520" class="mean"/><text x="388" y="524" class="tick">mean</text>',
        f'<line x1="445" y1="520" x2="472" y2="520" class="median"/><text x="478" y="524" class="tick">median</text>',
        '<text x="820" y="524" class="tick">Bars show percent of all pairs</text>',
        '</svg>',
    ])
    (ROOT / "jaccard_distribution.svg").write_text("\n".join(parts) + "\n", encoding="utf-8")


def export_run(con: duckdb.DuckDBPyConnection, run_id: str) -> dict[str, Any]:
    stats = make_pairs(con, run_id)
    write_csv(con,
        "SELECT study_id,title,species,sample_count,human_sample_count,refmet_row_count,unique_refmet_name_count "
        "FROM study WHERE run_id=? AND qualifies ORDER BY study_id",
        [run_id], ROOT/"eligible_studies.csv",
        ["study_id","title","species","sample_count","human_sample_count","refmet_rows","unique_refmet_names"])
    write_csv(con,
        "SELECT study_id,refmet_name FROM study_refmet_name WHERE run_id=? ORDER BY study_id,name_key",
        [run_id], ROOT/"refmet_name_sets.csv", ["study_id","refmet_name"])
    write_csv(con,
        "SELECT study1_id,study2_id,intersection_count,union_count,jaccard FROM pairwise_jaccard "
        "WHERE run_id=? ORDER BY jaccard DESC,study1_id,study2_id",
        [run_id], ROOT/"pairwise_jaccard.csv",
        ["study1_id","study2_id","intersection_count","union_count","jaccard"])
    write_csv(con,
        "SELECT human_sample_study_count,studies_with_refmet_names,eligible_study_count,pair_count,mean_jaccard,"
        "population_variance,sample_variance,zero_overlap_pairs,median_jaccard,minimum_jaccard,maximum_jaccard "
        "FROM summary WHERE run_id=?",
        [run_id], ROOT/"summary.csv",
        ["human_sample_study_count","studies_with_refmet_names","eligible_study_count","pair_count",
         "mean_jaccard","population_variance","sample_variance","zero_overlap_pairs","median_jaccard","minimum_jaccard","maximum_jaccard"])
    make_plot(con, run_id, stats)
    return stats


def run_live(database: Path) -> str:
    started = datetime.now(timezone.utc)
    run_id = started.strftime(RUN_ID_FORMAT)
    con = duckdb.connect(str(database))
    wb = Workbench(timeout=180, cache_dir=None)
    errors: list[dict[str, str]] = []
    manifest: dict[str, Any] = {
        "run_id": run_id,
        "started_at": started.isoformat(),
        "source": "Metabolomics Workbench only",
        "metabo_search_version": "0.1.0 (vendored Workbench client)",
        "human_sample_rule": "Explicit human sample taxon or sample source; when sample taxon is absent, a single-species Homo sapiens summary is a fallback; explicit non-human taxonomy overrides the summary.",
        "refmet_rule": "Only non-empty Workbench refmet_name values; no fallback to metabolite_name.",
        "name_normalization": "Unicode NFKC; collapse whitespace; casefold keys; one unique name per study/key.",
        "threshold": f"At least {MIN_REFMET_NAMES} unique RefMet names per study.",
        "variance": "Population variance (ddof=0) across all unique unordered pairs; sample variance also saved.",
        "workers": WORKERS,
    }
    con.execute(SCHEMA)
    con.execute("INSERT INTO run_info(run_id,started_at,status,metabo_search_version,errors_json,manifest_json) VALUES (?,?,?,?,?,?)",
                [run_id,started.isoformat(),"running",manifest["metabo_search_version"],"[]",json.dumps(manifest)])
    try:
        summaries = wb.refresh_corpus()
        all_summaries = list(summaries.values())
        candidates = [s for s in all_summaries if human_label(s.species)]
        manifest.update({"workbench_summary_count":len(all_summaries),"human_candidate_count":len(candidates)})
        con.execute("UPDATE run_info SET workbench_summary_count=?,human_candidate_count=?,manifest_json=? WHERE run_id=?",
                    [len(all_summaries),len(candidates),json.dumps(manifest),run_id])
        studies_by_id = {s.study_id:s for s in candidates}
        base_rows = [(run_id,s.study_id,s.study_title,s.species,"pending",None,None,"","","pending",None,None,False,"") for s in candidates]
        to_arrow_insert(con,"study",("run_id","study_id","title","species","sample_status","sample_count","human_sample_count","human_detection_method","sample_error","metabolite_status","refmet_row_count","unique_refmet_name_count","qualifies","metabolite_error"),base_rows)

        human_ids: list[str] = []
        for offset in range(0,len(candidates),BATCH_SIZE):
            chunk=candidates[offset:offset+BATCH_SIZE]
            details=list(wb.inspect(chunk,workers=WORKERS,include_summary=False,include_analyses=False,
                                    include_metabolites=False,on_error="collect"))
            failed=[d.study_id for d in details if d.error]
            if failed:
                retry=retry_failed_samples(wb,failed)
                retry_map={d.study_id:d for d in retry}
                details=[retry_map.get(d.study_id,d) if d.error else d for d in details]
            sample_records=[]
            for detail in details:
                summary=studies_by_id[detail.study_id]
                humans,method=find_human_samples(detail,summary)
                sample_status="error" if detail.error else "success"
                if detail.error:
                    errors.append({"study_id":detail.study_id,"phase":"samples","error":detail.error})
                    con.execute("UPDATE study SET sample_status='error',sample_error=?,human_detection_method='unknown' WHERE run_id=? AND study_id=?",
                                [detail.error,run_id,detail.study_id])
                    continue
                con.execute("UPDATE study SET sample_status='success',sample_count=?,human_sample_count=?,human_detection_method=?,metabolite_status=? WHERE run_id=? AND study_id=?",
                            [len(detail.samples),len(humans),method,"pending" if humans else "not_human",run_id,detail.study_id])
                if humans:
                    human_ids.append(detail.study_id)
                    sample_records.extend((run_id,detail.study_id,ordinal,name,organism,source,evidence,detection)
                                          for ordinal,name,organism,source,evidence,detection in humans)
            to_arrow_insert(con,"human_sample",("run_id","study_id","sample_ordinal","sample_name","organism_evidence","sample_source","evidence","detection_method"),sample_records)
            print(f"Sample verification: {min(offset+len(chunk),len(candidates))}/{len(candidates)}",flush=True)

        all_name_rows=[]
        for offset in range(0,len(human_ids),BATCH_SIZE):
            chunk=human_ids[offset:offset+BATCH_SIZE]
            details=list(wb.inspect(chunk,workers=WORKERS,include_summary=False,include_samples=False,
                                    include_analyses=False,include_metabolites=True,on_error="collect"))
            failed=[d.study_id for d in details if d.error]
            if failed:
                retry=retry_failed_metabolites(wb,failed)
                retry_map={d.study_id:d for d in retry}
                details=[retry_map.get(d.study_id,d) if d.error else d for d in details]
            for detail in details:
                if detail.error:
                    errors.append({"study_id":detail.study_id,"phase":"metabolites","error":detail.error})
                    con.execute("UPDATE study SET metabolite_status='error',metabolite_error=? WHERE run_id=? AND study_id=?",
                                [detail.error,run_id,detail.study_id])
                    continue
                names:dict[str,str]={}
                refmet_rows=0
                for row in detail.metabolites:
                    clean=canonical_name(row.refmet_name)
                    if clean is None:
                        continue
                    refmet_rows+=1
                    old=names.get(clean[0])
                    if old is None or clean[1]<old:
                        names[clean[0]]=clean[1]
                qualifies=len(names)>=MIN_REFMET_NAMES
                con.execute("UPDATE study SET metabolite_status='success',refmet_row_count=?,unique_refmet_name_count=?,qualifies=? WHERE run_id=? AND study_id=?",
                            [refmet_rows,len(names),qualifies,run_id,detail.study_id])
                if qualifies:
                    all_name_rows.extend((run_id,detail.study_id,key,name) for key,name in names.items())
            print(f"RefMet extraction: {min(offset+len(chunk),len(human_ids))}/{len(human_ids)}",flush=True)
        to_arrow_insert(con,"study_refmet_name",("run_id","study_id","name_key","refmet_name"),all_name_rows)
        stats=export_run(con,run_id)
        completed=datetime.now(timezone.utc).isoformat()
        manifest.update({"completed_at":completed,"human_sample_study_count":len(human_ids),
                         "eligible_study_count":stats["eligible_study_count"],"pair_count":stats["pair_count"],
                         "mean_jaccard":stats["mean_jaccard"],"population_variance":stats["population_variance"],
                         "sample_variance":stats["sample_variance"],"errors":errors,
                         "http_stats":wb.stats.__dict__})
        con.execute("UPDATE run_info SET completed_at=?,status='complete',human_sample_study_count=?,eligible_study_count=?,pair_count=?,mean_jaccard=?,population_variance=?,sample_variance=?,errors_json=?,manifest_json=? WHERE run_id=?",
                    [completed,len(human_ids),stats["eligible_study_count"],stats["pair_count"],stats["mean_jaccard"],
                     stats["population_variance"],stats["sample_variance"],json.dumps(errors),
                     json.dumps(manifest,ensure_ascii=False),run_id])
        (ROOT/"manifest.json").write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
        con.execute("CHECKPOINT")
        mean_text = f"{stats['mean_jaccard']:.10f}" if stats['mean_jaccard'] is not None else "n/a"
        variance_text = f"{stats['population_variance']:.10f}" if stats['population_variance'] is not None else "n/a"
        print(f"Complete: {stats['eligible_study_count']} eligible studies, {stats['pair_count']:,} pairs; "
              f"mean={mean_text}, population variance={variance_text}",flush=True)
        return run_id
    except Exception as exc:
        errors.append({"phase":"fatal","error":f"{type(exc).__name__}: {exc}"})
        con.execute("UPDATE run_info SET status='failed',completed_at=?,errors_json=?,manifest_json=? WHERE run_id=?",
                    [datetime.now(timezone.utc).isoformat(),json.dumps(errors),json.dumps(manifest),run_id])
        raise
    finally:
        wb.close()
        con.close()


def recompute(database: Path, run_id: str) -> None:
    con=duckdb.connect(str(database))
    try:
        stats=export_run(con,run_id)
        con.execute("UPDATE run_info SET eligible_study_count=?,pair_count=?,mean_jaccard=?,population_variance=?,sample_variance=?,status='complete' WHERE run_id=?",
                    [stats["eligible_study_count"],stats["pair_count"],stats["mean_jaccard"],
                     stats["population_variance"],stats["sample_variance"],run_id])
        row=con.execute("SELECT manifest_json FROM run_info WHERE run_id=?",[run_id]).fetchone()
        manifest=json.loads(row[0]) if row and row[0] else {"run_id":run_id}
        manifest.update(stats)
        manifest["recomputed_from_saved_database"] = True
        (ROOT/"manifest.json").write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
        con.execute("CHECKPOINT")
        mean_text = f"{stats['mean_jaccard']:.10f}" if stats['mean_jaccard'] is not None else "n/a"
        variance_text = f"{stats['population_variance']:.10f}" if stats['population_variance'] is not None else "n/a"
        print(f"Recomputed {stats['pair_count']:,} pairs; mean={mean_text}; "
              f"population variance={variance_text}")
    finally:
        con.close()


def main() -> None:
    parser=argparse.ArgumentParser(description="Workbench-only human study Jaccard analysis using RefMet names only.")
    parser.add_argument("--database",type=Path,default=DB_DEFAULT,help="DuckDB output (default: workbench_refmet_only.duckdb)")
    parser.add_argument("--recompute",metavar="RUN_ID",help="Recompute scores and outputs from a saved snapshot, without API calls")
    args=parser.parse_args()
    if args.recompute:
        recompute(args.database,args.recompute)
    else:
        run_live(args.database)


if __name__ == "__main__":
    main()
