# Metabo-search research

Reproducible findings from public metadata in MetaboLights and Metabolomics Workbench. The two analyses answer different questions: a cross-repository census of rows in human-tagged studies, and a Workbench-only comparison of RefMet-name sets.

## Finding 1: rows in human-tagged studies

In the checked-in census snapshot, **2,911 studies tagged as human** contain **517,895 parsed sample/factor rows**. Adding study-index counts for studies whose sample metadata could not be read gives an approximate total of **518,744 rows**. A broad source/tissue-text rule matched blood, plasma, or serum for **290,489 parsed rows** (approximately **290,989** with the same limited missing-metadata adjustment, about **56%** of the estimated row total).

| Repository | Studies tagged human | Parsed sample/factor rows | Rows matching blood-label rule | Estimated rows | Estimated blood-label rows |
|---|---:|---:|---:|---:|---:|
| MetaboLights | 1,038 | 214,334 | 106,714 | 215,138 | 107,214 |
| Metabolomics Workbench | 1,873 | 303,561 | 183,775 | 303,606 | 183,775 |
| **Combined** | **2,911** | **517,895** | **290,489** | **518,744** | **290,989** |

These are repository rows associated with studies tagged as human—not unique people, and not a sample-by-sample validation that every row is a human biological specimen. The saved census JSON has no acquisition timestamp, and repository indexes change over time. See [`reproductions/human-blood-stats/`](reproductions/human-blood-stats/README.md) for methods and limitations. In particular, the missing-metadata adjustments are study-level approximations; the blood-label rule is a broad text match.

## Finding 2: Workbench RefMet-name overlap

The Workbench-only analysis compares normalized distinct `refmet_name` sets across **734 human-evidenced studies with at least 100 such names**. Among the **269,011** unique unordered study pairs, the mean Jaccard overlap is **8.53%** and the median is **5.99%**; **13.47%** of pairs share no names.

Jaccard is the size of the intersection divided by the size of the union of two studies' recorded name sets. These figures describe name-set overlap in this selected, high-annotation subset. They do not estimate metabolite abundance or assert analytical or biological comparability. See [`reproductions/workbench-refmet-overlap/`](reproductions/workbench-refmet-overlap/README.md) for selection criteria, the saved snapshot, and caveats.

## Reproduce from scratch

Python 3.11 or newer is required. From the repository root, install the census dependencies and run the live client acquisition:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python reproductions/human-blood-stats/reproduce.py \
  --workers 4 --out reproductions/human-blood-stats/out
```

This queries the public APIs and typically takes around 13 minutes. To acquire the RefMet finding from the live Workbench API as well, install its pinned dependencies and run its acquisition script:

```bash
python -m pip install -r reproductions/workbench-refmet-overlap/requirements.txt
python reproductions/workbench-refmet-overlap/run_refmet_only.py
```

The RefMet script creates a fresh run in its DuckDB database and refreshes the CSV exports, manifest, and plot. Both live findings may change as upstream records change.

## Offline verification

The offline smoke tests validate the saved census outputs and RefMet DuckDB snapshot; they make no API requests:

```bash
python -m unittest discover -s tests -v
```

To recompute the RefMet pair scores and exports from the included snapshot without making API calls, use the run ID in its `manifest.json`:

```bash
python reproductions/workbench-refmet-overlap/run_refmet_only.py \
  --recompute 20261021T111615Z
```

## Provenance and license

Research code is provided under the MIT License. The census depends on the separately maintained [Metabo-search](https://github.com/bmds-tue/metabo-search) project at the immutable commit pinned in `requirements.txt`. The RefMet reproduction vendors a small client subset and preserves its upstream license in its `vendor/` directory. Repository data and metadata remain attributable to the respective public repositories.

## Publish to GitHub

Create an empty GitHub repository named `metabo-search-research`, then add it as this repository's remote and push the `main` branch:

```bash
git remote add origin https://github.com/<YOUR-ACCOUNT>/metabo-search-research.git
git push -u origin main
```
