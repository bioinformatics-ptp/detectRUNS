---
name: cran-check
description: Prepare and verify the detectRUNS R package for CRAN submission. Checks version (no dev suffix), date (today), NEWS.md entry, cran-comments.md, then runs R CMD check --as-cran. Use when the user wants to submit to CRAN, prepare a release, or run pre-CRAN checks.
tools: Read, Edit, PowerShell, Glob, Bash
---

# CRAN Submission Checker

Verify and fix the detectRUNS package for CRAN submission, then run `R CMD check --as-cran`.

The R package source is in the `detectRUNS/` subdirectory of the repo root.

## Step 1 — Get today's date

Use R to get today's date in YYYY-MM-DD format (cross-platform):

```r
Rscript -e "cat(format(Sys.Date(), '%Y-%m-%d'), '\n')"
```

## Step 2 — Check DESCRIPTION

Read `detectRUNS/DESCRIPTION` and verify:

### Version

The `Version:` field must NOT have a development suffix. Development versions look like:
- `0.9.6.9000`, `0.9.6.9001`, `0.9.6.9005`, etc. (fourth component ≥ 9000)

If the version has a dev suffix, remove it so the result is `<major>.<minor>.<patch>` only.

Example: `0.9.6.9005` → `0.9.7` (increment patch when stripping dev suffix, per r-pkgs.org convention).

Ask the user to confirm the target version before editing if unclear.

### Date

The `Date:` field must equal today's date in `YYYY-MM-DD` format. Update it if it differs.

## Step 3 — Check NEWS.md

Read `detectRUNS/NEWS.md`. There must be a top-level entry for the release version determined in Step 2, in the format:

```
# detectRUNS <version>
```

If the entry is missing, report it to the user and ask them to add the relevant changelog entries before proceeding. Do NOT auto-generate NEWS content — changelog entries describe user-visible changes and require human judgement.

## Step 4 — Check cran-comments.md

Read `detectRUNS/cran-comments.md`. Remind the user to update it with current R CMD check results after the check completes (Step 5). It should include:
- Test environments (OS / R version matrix)
- R CMD check result summary (`N errors | N warnings | N notes`)
- Explanation for any NOTEs
- Downstream dependencies section

## Step 5 — Run R CMD check --as-cran

Run the check from the repo root using devtools, which sets the correct CRAN-equivalent flags (`--as-cran` enables `NOTE`-level checks for URLs, examples, etc.):

```sh
Rscript -e "devtools::check('detectRUNS', remote = TRUE, manual = TRUE, error_on = 'never')"
```

- `remote = TRUE` — checks for packages only available on CRAN (not just locally)
- `manual = TRUE` — builds the PDF manual (requires LaTeX; skip with `manual = FALSE` if LaTeX is not available)
- `error_on = 'never'` — show full output even if there are errors/warnings

Capture and display the full output. Summarise:
- Number of ERRORs, WARNINGs, NOTEs
- Text of each ERROR/WARNING/NOTE

## Step 6 — Final report

After all steps, produce a checklist summary:

```
CRAN Submission Readiness
─────────────────────────
[ ] Version: <value> — no dev suffix
[ ] Date:    <value> — matches today
[ ] NEWS.md: entry for <version> present
[ ] R CMD check: N errors | N warnings | N notes
[ ] cran-comments.md: updated
```

Mark each item ✅ (OK) or ❌ (needs attention). List any blocking items (errors, warnings, missing NEWS entry) that must be resolved before submission.

## Important notes

- **ERRORs and WARNINGs are blocking** — CRAN will reject the package.
- **NOTEs should be minimised** — if unavoidable (e.g. installed size > 5Mb due to compiled code), document them in `cran-comments.md`.
- The known NOTE about installed package size (7.1Mb due to Rcpp + example data) is pre-documented in `cran-comments.md` and is acceptable.
- After a successful CRAN submission and acceptance, bump the version to the next dev version with `usethis::use_dev_version()`.
