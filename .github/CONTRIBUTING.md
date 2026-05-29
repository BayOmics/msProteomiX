# Contributing to msProteomiX

## Development Setup

```bash
git clone https://github.com/BayOmics/msProteomiX.git
cd msProteomiX

# Install dependencies
Rscript -e 'devtools::install_deps(".", dependencies = TRUE)'

# Load for development
Rscript -e 'devtools::load_all(".")'
```

## Code Standards

### File I/O
- **Always use `base::read.delim()`** for TSV/CSV files
- **Never use `readr::read_tsv()` or `vroom`** — they cause memory corruption on Apple Silicon RStudio
- Use `.read_omics_file()` from `parsers.R` for all data file reading

### Encoding
- **Runtime strings must be ASCII only** — no emoji, no em-dash, no CJK characters in `message()`, `stop()`, `cat()`, `warning()`, `paste()`
- **Comments may contain CJK** — R ignores comments at runtime
- Use Unicode escapes (`\u274c`) only if absolutely needed, prefer English text

### NAMESPACE
- All exported functions need `#' @export` roxygen tag
- After adding exports: `devtools::document()`
- NSE column names → add to `R/globals.R` via `utils::globalVariables()`
- New package imports → add to `R/imports.R` via `#' @importFrom`

### Package Dependencies
- **Imports**: Core packages always available (dplyr, ggplot2, etc.)
- **Suggests**: Optional packages checked with `requireNamespace()` before use (limma, Biostrings, etc.)

## Adding a New Search Engine

See [design_specification.md](../vignettes/design_specification.md) Section 3.2 for the template.

## Commit Messages

Follow [Conventional Commits](https://www.conventionalcommits.org/):

```
feat: add new feature
fix: fix a bug
refactor: code restructuring (no behavior change)
docs: documentation only
chore: maintenance (gitignore, CI, etc.)
perf: performance improvement
```

## R CMD Check

Before submitting changes, ensure:

```bash
Rscript -e 'devtools::check(".", args="--no-examples --no-tests --no-vignettes")'
# Target: 0 errors, 0 warnings
```

## Known Pitfalls

> These are the most common mistakes when developing msProteomiX.
> See `vignettes/design_specification.md` Section 5 for the full list.

1. **`.msProteomiX_env` must be `@export`ed** — user scripts reference it directly
2. **`topTable(sort.by = "none")`** — Limma default sort misaligns rows with input matrix
3. **`stat_ellipse` requires ≥3 samples** — always guard with `if (n >= 3)`
4. **FragPipe quantification columns** — regex must anchor `MaxLFQ Intensity$` (suffix, not prefix)
5. **Spectronaut tibble** — `new_MsDataSet()` forces `as.data.frame()` but verify this persists
6. **`devtools::install()` while package is loaded** — corrupts lazy-load DB. Always restart R first.
