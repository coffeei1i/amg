# Server quick start

Copy this directory to any server location, then run from its parent directory.
No biological input path is embedded in the public code.

```bash
cd BRS-Map

Rscript scripts/run_brs_map.R --stage preflight
```

For a formal stage, copy `examples/server_adapter_template.R` to a private
working directory, implement only its data-loading and checkpoint-writing
functions, and run:

```bash
Rscript scripts/run_brs_map.R \
  --stage all \
  --input-root /your/input/root \
  --output-root /your/new/output/root \
  --adapter /your/private/brs_map_adapter.R
```

Use `--resume` only when the adapter performs strict checkpoint validation.
Without `--resume`, BRS-Map refuses to write into a non-empty formal output
directory.

The public release deliberately separates stable mathematics from private I/O:
edit the adapter for server object paths; do not duplicate or alter the frozen
thresholds in `R/contracts.R`.
