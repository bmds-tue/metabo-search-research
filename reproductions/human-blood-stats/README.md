# Repository rows in human-tagged metabolomics studies
*Client-only reproduction (no database).*

**Finding.** In the checked-in metadata snapshot, **2,911 studies tagged as human** contain **517,895 parsed sample/factor rows**. Adding study-index counts for studies whose sample metadata could not be read gives an approximate total of **518,744 rows**. A broad tissue/source-text rule matched blood, plasma, or serum for **290,489 parsed rows** (approximately **290,989** with the same limited missing-metadata adjustment; about **56%** of the estimated row total).

These are repository rows associated with human-tagged studies—not unique people, and not a sample-by-sample validation that every row is a human biological specimen. Counts include whatever rows the repositories expose, such as repeated measurements, technical replicates, blanks, and QC material.

## Saved snapshot results

| repository | studies tagged human | parsed sample/factor rows | rows matching blood-label rule | estimated rows | estimated blood-label rows |
|---|---:|---:|---:|---:|---:|
| MetaboLights | 1,038 | 214,334 | 106,714 | 215,138 | 107,214 |
| Metabolomics Workbench | 1,873 | 303,561 | 183,775 | 303,606 | 183,775 |
| **Combined** | **2,911** | **517,895** | **290,489** | **518,744** | **290,989** |

The checked-in census JSON does not record its acquisition timestamp; repository indexes can change. Treat these as counts from the saved snapshot, not live guarantees. The estimate adds study-index sample counts for 8 MetaboLights studies (804 rows) and 2 Workbench studies (45 rows) with missing/empty sample-level metadata. For MetaboLights, the blood-label estimate adds all rows from an unreadable study when a study-level facet contains a matching term. For Workbench, no blood labels are inferred for rows missing factor metadata, so its blood-label total remains the parsed count.

## Method

The two public repositories are queried independently; there is no cross-repository join or deduplication.

- **MetaboLights.** `search(organisms=["Homo sapiens"])` selects studies tagged with the human organism. Each selected study is inspected for its ISA sample sheet. The script counts every parsed row and applies the blood-label rule to `organism_part`.
- **Workbench.** The whole-index summary corpus is filtered to `species == "Homo sapiens"`. The script fetches each selected study's factor rows, counts every row, and applies the rule to `sample_source`.

Thus “human” is the repository's study-level organism/species annotation; the census does not reclassify each sample row by species. The blood-label rule is case-insensitive and matches text containing `blood`, `plasma`, or `serum`, except terms containing `seminal`.

## Interpretation and limitations

- The unit is a repository sample/factor row, not a unique participant. Repeated measures and technical replicates count separately; blank and QC rows may also count.
- Study-level human tagging does not guarantee that every row in a study is a human biospecimen. Cell lines, controls, or other row types remain in the count.
- The blood figure is a broad text-label tally, not a curated tissue classification. Different vocabularies and fields are used by the two repositories; ambiguous labels may be included or missed.
- The missing-metadata adjustments are approximate study-level estimates, not individually classified samples. In particular, the MetaboLights blood fallback assigns the full study index count to blood when any study-level facet matches.
- The repositories are counted separately, with no cross-repository deduplication.

## Reproducibility

Install the pinned dependencies from the repository root, then run from that root:

```bash
python -m pip install -r requirements.txt
python reproductions/human-blood-stats/reproduce.py --workers 4
```

- **Library:** `metabo-search` at commit `4ee8d01525a898b83dfcd3ec408e372f013236f8`.
- **Runtime:** about 13 minutes; the run is network-bound.
- **Outputs:** `out/metabolights.json`, `out/workbench.json`, and `out/summary.json` (counts plus per-study detail for auditing).
