#!/usr/bin/env bash
# R package coverage for rsxr via covr (R/ + src/rsxr.c when instrumentable).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
OUT="${1:-r_coverage.xml}"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript not found" >&2
  exit 1
fi

echo "==> stage rsxcore into rsx-r"
sh "$ROOT/rsx-r/tools/stage_rsxcore.sh" "$ROOT/rsx-r"

echo "==> covr::package_coverage(rsx-r)"
Rscript -e "
options(warn = 1)
if (!requireNamespace('covr', quietly = TRUE)) {
  install.packages('covr', repos = 'https://cloud.r-project.org')
}
# Ensure package deps for tests
for (p in c('testthat', 'tibble', 'readr')) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = 'https://cloud.r-project.org')
  }
}
cov <- covr::package_coverage(
  path = 'rsx-r',
  type = 'tests',
  quiet = FALSE,
  clean = TRUE
)
print(cov)
covr::to_cobertura(cov, filename = '${OUT}')
cat('OK wrote ${OUT}\\n')
"
test -s "$OUT"
echo "OK R coverage $OUT"
