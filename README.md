# mheatmap

Experimental collaboration repository for preparing the VIS 2026 short paper for `Q9gJYx/mheatmap`.

## Repository structure

- `src/`: source code and reusable modules
- `data/`: dataset-related files (keep large/private data out of Git)
- `output/`: generated artifacts, figures, and intermediate results
- `docs/`: project notes, paper-prep documentation, and writing support files

## Environment management (uv)

1. Install [uv](https://docs.astral.sh/uv/)
2. Sync the project environment:

   ```bash
   uv sync
   ```

3. Run commands in the managed environment:

   ```bash
   uv run python -V
   ```

## Notes

- Commit `uv.lock` when dependencies are introduced to keep environments reproducible.
- Keep generated outputs and large data files outside source control unless explicitly needed.
