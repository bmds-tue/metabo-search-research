# Human and blood-derived sample counts in MetaboLights and Metabolomics Workbench
*Client-only reproduction (no database).*

**Finding.** Across the two major public metabolomics repositories, the
MetaboLights and Metabolomics Workbench studies with a human organism contain
**≈ 512,000 samples**, of which **≈ 285,000 (~56%) are blood-derived**
(blood, plasma, or serum). This is a per-sample census of the repositories'
own metadata, not a curated subset.

This reproduction uses only the library clients: it fetches each repository's
sample metadata and counts rows directly in Python.

## Data sources

| repository | studies (all) | human studies | human samples | blood-derived |
|---|---:|---:|---:|---:|
| MetaboLights | 3,511 | 1,025 | 213,952 | 106,704 |
| Metabolomics Workbench | 4,629 | 1,871 | 297,536 | 178,201 |
| **combined** | **8,140** | **2,896** | **511,488** | **284,905** |

MetaboLights samples are read from the ISA-Tab sample sheets of each human
study; Workbench samples are read from each study's `factors` endpoint. A
sample is counted as blood-derived when its tissue/source contains
`blood`, `plasma`, or `serum`; `seminal plasma` is excluded. Counts are
per-sample rows, not unique individuals.

## Method

No database and no working set — the two clients are called directly and the
rows are counted in Python.

- **MetaboLights.** `search(organisms=["Homo sapiens"])` retrieves the human
  study set; each study is deep-inspected for its ISA sample sheet
  (`inspect(include_metabolites=False)`), and samples are classified on
  `organism_part`.
- **Workbench.** The whole-index summary corpus is filtered to
  `species == "Homo sapiens"`; each study's factor rows are fetched
  (`samples(study_id)`) and classified on `sample_source`.

Both repositories are counted independently; no cross-repository join is
performed.

## Results

**≈ 512,000 human samples; ≈ 285,000 (~56%) blood-derived.**

Blood is dominated by the canonical biofluids: in MetaboLights, `blood plasma`
(45,790) and `blood serum` (33,366) account for ~74% of the blood-derived
total, with the remainder spread over `plasma`/`serum`/`blood` variants, whole
and cord blood, and blood cells. In the Workbench, `Blood (plasma)` (70,561)
and `Blood (serum)` (28,837) similarly dominate, followed by lowercase
`plasma` (24,690) and `Blood` (13,845).

For full coverage, the estimate adds the search-index `sample_count` for the
studies whose per-sample metadata could not be read (8 MetaboLights studies
without ISA files on the mirror; 2 Workbench studies whose factor endpoint is
empty), giving **512,337 samples / 285,405 blood**.

## Limitations

- **Sample rows, not individuals.** Repeated sampling of the same subject
  (time points, technical replicates) is counted as separate samples; QC,
  blank, and solvent rows are included in the totals.
- **Live repositories.** Both indexes change continuously, so the counts are a
  dated snapshot and move by a few studies between runs.
- **Unreadable metadata.** 10 studies (~0.3% of the corpus) expose no
  per-sample metadata and are covered only by their study-level sample count.
- **Vocabulary.** Blood classification is substring-based on the repositories'
  free-text tissue/source fields; borderline terms (e.g. bone-marrow serum,
  dried blood spots) are included under a broad "blood of some sort" reading.

## Reproducibility

Install the pinned dependencies once from the repository root, then run this
script from the repository root:

```bash
python -m pip install -r requirements.txt
python reproductions/human-blood-stats/reproduce.py --workers 4
```

- **Library:** `metabo-search` at commit
  `4ee8d01525a898b83dfcd3ec408e372f013236f8`.
- **Runtime:** ~13 minutes; the run is network-bound. `--workers 4` is the
  measured optimum (8 triggers server throttling, 2 underutilizes).
- **Output:** `reproductions/human-blood-stats/out/metabolights.json`,
  `out/workbench.json`, and `out/summary.json` (counts plus per-study detail
  for auditing).
