# Contributing to clustOpt

## Before You Start

- **Open an issue first.** Before writing code, open an issue describing the bug or feature so we can discuss the approach. This avoids wasted effort on changes that may not fit the project direction.
- **Keep PRs focused.** One bug fix or feature per PR. Smaller PRs are easier to review and merge.
- **Run `devtools::check()` before submitting.** CI can be slow with the R dependency install -- catching issues locally saves time for everyone.
- **Questions?** Open a [GitHub issue](https://github.com/gladstone-institutes/clustOpt/issues) if you get stuck or need guidance.

## Development Workflow

1. [Fork the repository and create a feature branch](https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project) from `main`.
2. Make code changes in `R/`, `man/`, etc.
3. `devtools::document()` to regenerate roxygen docs.
4. `devtools::test()` to run tests.
5. When ready to merge, rebuild the site (includes the vignette):
   ```r
   pkgdown::build_site()
   ```
6. Commit everything including `docs/`.
7. Open a pull request against `main`. CI will run R-CMD-check and verify the vignette compile date is current.

## Development Setup

```r
# Install dependencies
devtools::install_deps(dependencies = TRUE)

# Load package for interactive development
devtools::load_all()

# Run tests
devtools::test()

# Check package
devtools::check()
```

## Code Style

This package follows the [tidyverse style guide](https://style.tidyverse.org/). The [`styler`](https://styler.r-lib.org/) package can help format code automatically:

```r
# Style a single file
styler::style_file("R/clust_opt.R")

# Style the whole package
styler::style_pkg()
```

Use `roxygen2` for documentation and run `devtools::document()` to regenerate man pages.

Verbosity levels for user-facing functions: 0 = silent, 1 = milestones, 2 = detailed, 3 = includes Seurat output. Use `cli` for all user-facing messages.

## Versioning

This package uses [semantic versioning](https://semver.org/) (`MAJOR.MINOR.PATCH`). Most PRs should **not** bump the version -- version bumps are coordinated by maintainers at release time.

- **PATCH** (e.g., 1.2 -> 1.2.1): Bug fixes, documentation updates, internal refactors with no user-facing changes.
- **MINOR** (e.g., 1.2 -> 1.3): New features, new exported functions, or changes to default behavior.
- **MAJOR** (e.g., 1.3 -> 2.0): Breaking changes to the API (renamed/removed functions, changed return types).

When a version bump is needed, update `Version` and `Date` in `DESCRIPTION` and add a section to `NEWS.md`.

## Rebuilding the Vignette and Site

The vignette requires the tutorial dataset and ~32GB RAM, so it cannot be built in CI. Before merging any PR that changes package source or documentation, rebuild the site locally:

```r
pkgdown::build_site()
```

This rebuilds the full site including reference docs and the vignette.

Commit the updated `docs/` directory with your PR.

The CI will warn on PRs if the vignette's compile date is older than the most recent source change. The check uses day granularity, so rebuild after all source changes are finalized -- any source commits made after the vignette was last rendered will trigger the warning.
