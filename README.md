# Metabo-search research

Reproducible research analyses built on public metadata from MetaboLights and
Metabolomics Workbench. The repository contains one client-based reproduction
of a human/blood sample census and a separate Workbench RefMet-overlap
analysis.

## Human and blood-derived sample counts

The census asks how many samples belong to human studies, and how many are
blood-derived (blood, plasma, or serum). It counts sample rows, not unique
people, and excludes `seminal plasma`.

| Repository | Human studies | Parsed sample rows | Blood-derived rows |
|---|---:|---:|---:|
| MetaboLights | 1,025 | 213,952 | 106,704 |
| Metabolomics Workbench | 1,871 | 297,536 | 178,201 |
| **Combined** | **2,896** | **511,488** | **284,905** |

Eight MetaboLights studies and two Workbench studies have no readable sample
rows. Adding their study-index sample counts gives the full-coverage estimate
of **512,337 samples**, of which **285,405 are estimated to be blood-derived**.
The estimates are approximate because the missing studies cannot be classified
at sample level. Repository indexes change over time; these are saved snapshot
results, not live guarantees.

See the [client-only reproduction](reproductions/human-blood-stats/README.md)
for the method, per-study audit outputs, and instructions to rerun it. It counts
rows directly from the two clients without materializing a full working-set
database.

## Additional analysis: Workbench RefMet overlap

[`analyses/workbench-refmet-only/`](analyses/workbench-refmet-only/README.md)
compares normalized RefMet-name sets for human-sample-verified Workbench
studies. Its saved 21 October 2026 DuckDB snapshot supports exact offline
recomputation of the pair scores, exports, and plot. This is a distinct
analysis and does not query MetaboLights.

## Setup and reproduction

Python 3.11 or newer is required. The client-only census uses Metabo-search at
the pinned commit listed in `requirements.txt`.

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python reproductions/human-blood-stats/reproduce.py --workers 4
```

The live run queries the public APIs and takes around 13 minutes. It writes
per-study JSON and summary outputs to the reproduction's `out/` directory.
Results can differ as upstream records change.

The auxiliary RefMet analysis has its own vendored client and dependencies; see
its README for its live-run and offline `--recompute` commands.

## Offline checks

The offline smoke suite validates the saved census results and classification
rules. It does not repeat the API acquisition:

```bash
python -m unittest discover -s tests -v
```

The RefMet analysis can be recomputed offline from its included database:

```bash
cd analyses/workbench-refmet-only
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python run_refmet_only.py --recompute 20261021T074807Z
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
analysis vendors a small client subset and preserves that upstream license in
its `vendor/` directory. Repository data and metadata remain attributable to
the respective public repositories.

## Publish to GitHub

Create an empty GitHub repository named `metabo-search-research`, then add it as
this repository's remote and push the `main` branch:

```bash
git remote add origin https://github.com/<YOUR-ACCOUNT>/metabo-search-research.git
git push -u origin main
```
