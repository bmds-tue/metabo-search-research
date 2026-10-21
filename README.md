# Metabo-search research

Reproducible findings from public metadata in MetaboLights and Metabolomics
Workbench. Both findings live under [`reproductions/`](reproductions/): a
human/blood sample census and a distinct Workbench RefMet-overlap analysis.

## Finding 1: human and blood-derived samples

The census asks how many samples belong to human studies, and how many are
blood-derived (blood, plasma, or serum). It counts sample rows, not unique
people, and excludes `seminal plasma`.

| Repository | Human studies | Parsed sample rows | Blood-derived rows |
|---|---:|---:|---:|
| MetaboLights | 1,038 | 214,334 | 106,714 |
| Metabolomics Workbench | 1,873 | 303,561 | 183,775 |
| **Combined** | **2,911** | **517,895** | **290,489** |

Eight MetaboLights studies and two Workbench studies have no readable sample
rows. Adding their study-index sample counts gives the full-coverage estimate
of **518,744 samples**, of which **290,989 are estimated to be blood-derived**.
The estimates are approximate because the missing studies cannot be classified
at sample level. Repository indexes change over time; these are saved snapshot
results, not live guarantees.

See [`reproductions/human-blood-stats/`](reproductions/human-blood-stats/README.md)
for the method, audit outputs, and rerun instructions. It counts directly from
the repository clients without materializing a full working-set database.

## Finding 2: Workbench RefMet-name overlap

[`reproductions/workbench-refmet-overlap/`](reproductions/workbench-refmet-overlap/README.md)
compares normalized RefMet-name sets for human-sample-verified Workbench
studies. The fresh snapshot contains 734 eligible studies and 269,011
unordered pairs, with mean Jaccard overlap **0.085345 (8.53%)**. Its DuckDB
snapshot supports exact offline recomputation of the pair scores, exports, and
plot. This is a distinct finding and does not query MetaboLights.

## Reproduce from scratch

Python 3.11 or newer is required. From the repository root, install the census
dependencies and run the live client acquisition:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python reproductions/human-blood-stats/reproduce.py \
  --workers 4 --out reproductions/human-blood-stats/out
```

This queries the public APIs and typically takes around 13 minutes. To acquire
the RefMet finding from the live Workbench API as well, install its pinned
dependencies and run its acquisition script:

```bash
python -m pip install -r reproductions/workbench-refmet-overlap/requirements.txt
python reproductions/workbench-refmet-overlap/run_refmet_only.py
```

The RefMet script performs a fresh live acquisition, adds a run to its DuckDB
database, and refreshes the CSV exports, manifest, and plot. Both findings may
change as upstream records change.

## Offline verification

The offline smoke tests validate the saved census outputs and the saved RefMet
DuckDB snapshot; they make no API requests:

```bash
python -m unittest discover -s tests -v
```

To recompute the RefMet pair scores and exports from the included snapshot
without making API calls, use the run ID in its `manifest.json`:

```bash
python reproductions/workbench-refmet-overlap/run_refmet_only.py \
  --recompute 20261021T111615Z
```

## Interpretation and limitations

- Counts are rows in repository metadata, not unique participants. Repeated
  measures, technical replicates, blanks, and QC samples can all contribute.
- Blood classification is a broad, case-insensitive substring rule applied to
  source/tissue text; borderline labels can be included.
- Unreadable studies are estimated from study-level sample counts and facets,
  not individual rows.
- The repositories are counted independently; there is no cross-repository
  deduplication or join.

## License and provenance

Research code is provided under the MIT License. The census depends on the
separately maintained [Metabo-search](https://github.com/bmds-tue/metabo-search)
project at the immutable commit pinned in `requirements.txt`. The RefMet
reproduction vendors a small client subset and preserves that upstream license
in its `vendor/` directory. Repository data and metadata remain attributable
to the respective public repositories.

## Publish to GitHub

Create an empty GitHub repository named `metabo-search-research`, then add it as
this repository's remote and push the `main` branch:

```bash
git remote add origin https://github.com/<YOUR-ACCOUNT>/metabo-search-research.git
git push -u origin main
```
