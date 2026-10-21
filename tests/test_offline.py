"""Offline smoke checks for the checked-in research snapshots."""
from __future__ import annotations

import importlib.util
import json
import unittest
from pathlib import Path

import duckdb

ROOT = Path(__file__).resolve().parents[1]


def load_script(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class HumanBloodSnapshotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.client = load_script(
            "client_reproduction",
            ROOT / "reproductions/human-blood-stats/reproduce.py",
        )
        cls.summary = json.loads(
            (ROOT / "reproductions/human-blood-stats/out/summary.json").read_text()
        )
        cls.workbench = json.loads(
            (ROOT / "reproductions/human-blood-stats/out/workbench.json").read_text()
        )

    def test_saved_parsed_counts(self) -> None:
        expected = {
            "metabolights": (1038, 214334, 106714),
            "workbench": (1873, 303561, 183775),
        }
        for repository, (studies, samples, blood) in expected.items():
            with self.subTest(repository=repository):
                result = self.summary[repository]
                self.assertEqual(result["human_studies"], studies)
                self.assertEqual(result["human_samples_total"], samples)
                self.assertEqual(result["human_samples_blood"], blood)

    def test_estimated_combined_totals_match_reported_finding(self) -> None:
        sample_estimate = sum(
            self.summary[repository]["human_samples_total_estimate"]
            for repository in ("metabolights", "workbench")
        )
        blood_estimate = sum(
            self.summary[repository]["human_samples_blood_estimate"]
            for repository in ("metabolights", "workbench")
        )
        self.assertEqual(sample_estimate, 518744)
        self.assertEqual(blood_estimate, 290989)

    def test_missing_workbench_rows_are_added_to_estimate(self) -> None:
        missing = sum(
            record["summary_samples"]
            for record in self.workbench["per_study"]
            if record["factor_rows"] == 0
        )
        self.assertEqual(missing, 45)
        self.assertEqual(
            self.workbench["human_samples_total_estimate"], 303561 + missing
        )
        self.assertEqual(
            self.summary["workbench"]["human_samples_total_estimate"], 303606
        )

    def test_blood_classification(self) -> None:
        self.assertTrue(self.client.is_blood("Blood (plasma)"))
        self.assertTrue(self.client.is_blood("blood cells"))
        self.assertFalse(self.client.is_blood("seminal plasma"))
        self.assertFalse(self.client.is_blood("urine"))
        self.assertTrue(self.client.is_canonical("blood plasma"))
        self.assertFalse(self.client.is_canonical("Blood (plasma)"))


class RefMetSnapshotTests(unittest.TestCase):
    def test_saved_duckdb_snapshot_matches_manifest_offline(self) -> None:
        analysis = ROOT / "reproductions/workbench-refmet-overlap"
        manifest = json.loads((analysis / "manifest.json").read_text())
        with duckdb.connect(str(analysis / "workbench_refmet_only.duckdb"), read_only=True) as db:
            row = db.execute(
                "SELECT eligible_study_count, pair_count, mean_jaccard, population_variance "
                "FROM summary WHERE run_id = ?",
                [manifest["run_id"]],
            ).fetchone()
        self.assertIsNotNone(row)
        self.assertEqual(row[0], 734)
        self.assertEqual(row[1], 269011)
        self.assertAlmostEqual(row[2], manifest["mean_jaccard"], places=14)
        self.assertAlmostEqual(row[3], manifest["population_variance"], places=14)


if __name__ == "__main__":
    unittest.main()
