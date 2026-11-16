## Rules

- Do NOT run `git commit`, `git push`, or any other git write commands. The human will handle all git operations.
- Do NOT create, read, or modify `.env` files or any file containing credentials, tokens, or secrets.
- Do NOT make network requests to external services (no `curl`, `wget`, API calls, etc.).
  - Exception: with the human's explicit approval for each request, you may download packages
    from package repositories (e.g., PyPI, CRAN, Bioconductor) into a throwaway environment.
    - A throwaway environment lives outside the repository (e.g., a `uv venv` under `/tmp`) and
      is used only for one-off inspection or prototyping.
    - Before asking, state the environment path, the packages, and why they're needed.
    - Approval covers only that request; ask again for any new packages or environments.
    - Do not use throwaway environments to run analysis code whose outputs are kept in the repo.
    - Throwaway environments must never be added to the project's renv or Python environment files.
  - Exception: with the human's explicit approval for each request, you may install packages
    into a module's managed environment: its renv library (`renv.lock`), its uv project
    (`pyproject.toml` + `uv.lock`), or a conda environment defined by an `environment.yml`
    in the module.
    - Before asking, state: the module, the environment, each package and version, its source
      (CRAN, Bioconductor, PyPI, a named conda channel such as conda-forge or bioconda, or
      GitHub with a commit), why it's needed, and which files will change (e.g.,
      `dependencies.R`, `renv.lock`, `pyproject.toml`, `uv.lock`, `environment.yml`).
    - GitHub or other non-repository sources need a stated reason why the package isn't
      available from CRAN, Bioconductor, PyPI, or conda-forge/bioconda. Do not add new conda
      channels without approval.
    - Only add what was approved. Do not upgrade, downgrade, or remove existing packages
      unless that change was approved explicitly; if an install would change other packages,
      stop and report which ones before continuing.
    - Install only through the module's tooling, never into the user or system library:
      - R: `renv::install()`, then `renv::snapshot()`.
      - Python (uv): `uv add`, which updates `pyproject.toml` and `uv.lock`; do not use
        `pip install` into the project environment.
      - Conda: edit `environment.yml`, then `conda env update --file environment.yml --prune`;
        do not `conda install` into an environment ad hoc or into `base`.
    - Afterwards, summarize the diff of the environment files, including any dependencies
      pulled in, and record the change in the relevant `.workspace/` plan.
- Do NOT install system-level packages or modify system configuration.
- Do NOT access or reference any data that is not already present in this repository.

## Expectations

- When you have made changes, stop and describe what you changed and why so the human can review.
- If you are unsure whether an action is allowed, ask rather than proceeding.
- Keep changes small and incremental. Prefer multiple small steps over one large batch of changes.

## Working on features with `.workspace/`

`.workspace/` holds planning notes and scratch documents for in-progress features. It is a place to think and plan, not a place for analysis code or data.

- Each feature gets its own subdirectory, named after the analysis module it belongs to (e.g., `.workspace/annotate-aml/` for work in `analyses/annotate-aml/`).
- Within that subdirectory, use one Markdown file per feature or task, with a descriptive name (e.g., `.workspace/annotate-aml/reference.md` for building the AML scAtlas reference).
- Before starting work on a feature, read its plan in `.workspace/` if one exists, and follow it.
- Write a plan to `.workspace/` before making substantial changes, and get the human's agreement on it first.
- Keep plans current: as steps are completed, decisions are made, or the approach changes, update the plan file so it reflects the actual state of the work. Record decisions and the reasons for them.
- Analysis code, notebooks, and outputs belong in `analyses/<module>/`, not in `.workspace/`.
