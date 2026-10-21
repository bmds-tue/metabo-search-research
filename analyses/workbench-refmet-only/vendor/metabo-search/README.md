# Vendored Workbench client subset

This contains only the modules needed to query and inspect Metabolomics Workbench using the typed `Workbench` client. It is copied from `bmds-tue/metabo-search` version 0.1.0, commit `176d70faa1e7690f5bae0c94dc9aac48389e4fc6` (MIT License). MetaboLights modules and unrelated ingestion helpers are intentionally omitted.

`run_refmet_only.py` adds this `src/` directory to `sys.path`, so no editable install or separate metabo-search installation is needed.
