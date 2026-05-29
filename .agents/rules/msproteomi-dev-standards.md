---
trigger: always_on
---

- Do NOT use readr, vroom, or data.table for file I/O. Use base R only (read.delim, read.csv).
- Do NOT edit NAMESPACE manually. Use devtools::document() to regenerate.
- Do NOT put non-ASCII characters in runtime strings (message/stop/cat/warning). Use \uXXXX escapes or English.
- Do NOT use Limma topTable() without sort.by = "none".
- Do NOT use `get(org_db)` to access Bioconductor annotation databases inside package functions. Use `loadNamespace(org_db)[[org_db]]`.
- Do NOT use `Label_Name` for enrichment gene extraction. `Label_Name` = UniProt Entry Name (e.g., "NUD4B_HUMAN"), not gene symbol. Always prioritize the "Gene" column.
- Wrap calls to optional Bioconductor packages (ReactomePA, STRINGdb) in tryCatch with a user-friendly install message, never let the script stop with an opaque error.
- Respond to user in Chinese. Commit messages in English (Conventional Commits).
- Always check the msProteomiX Knowledge Item before making changes.
- After modifying R source code, MUST run `devtools::install()` to update the installed package. `devtools::load_all()` only works in the dev session — user's RStudio loads the installed version.
- Keep .Rbuildignore up to date. Demo data, raw MS files, and large test datasets MUST be excluded. The package source dir and user's working dir (wkdir/) are independent — search engine output files must NEVER be bundled into the R package.
- When adding new user scripts to inst/scripts/, remember that `create_project()` copies from the **installed** package, not from the source tree. The user will only see new scripts after `devtools::install()`.