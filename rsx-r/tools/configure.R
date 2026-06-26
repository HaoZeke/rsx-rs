# Check that cargo is available and (optionally) new enough. The minimum
# supported Rust version is read from Config/rsxr/MSRV in DESCRIPTION.

safe_system2 <- function(cmd, args) {
  out <- tempfile()
  on.exit(unlink(out, force = TRUE))
  ret <- suppressWarnings(system2(cmd, args, stdout = out, stderr = out))
  list(success = identical(ret, 0L),
       output = if (file.exists(out)) readLines(out, warn = FALSE) else "")
}

get_msrv <- function() {
  x <- tryCatch(read.dcf("DESCRIPTION", fields = "Config/rsxr/MSRV")[[1]],
                error = function(e) NA_character_)
  if (length(x) == 0) NA_character_ else x
}

cat("*** Checking if cargo is installed\n")
res <- safe_system2("cargo", "version")
if (!isTRUE(res$success)) {
  cat("
-------------- ERROR: CONFIGURATION FAILED --------------------

The 'cargo' command is not available. rsxr compiles the rsxcore Rust
library, so a Rust toolchain is required.

Install Rust from <https://rustup.rs/> (or your distribution's packages)
and ensure 'cargo' and 'rustc' are on PATH.

The configure script also stages rsxcore into src/rust/rsxcore from
RSX_CORE_DIR or the monorepo sibling ../rsxcore when that tree is absent.

---------------------------------------------------------------

")
  quit("no", status = 2)
}

msrv <- get_msrv()
if (!is.na(msrv)) {
  m <- regmatches(res$output, regexec("cargo\\s+(\\d+\\.\\d+\\.\\d+)", res$output))[[1]]
  if (length(m) == 2 && package_version(m[2]) < package_version(msrv)) {
    cat(sprintf("ERROR: cargo %s is older than the required %s\n", m[2], msrv))
    quit("no", status = 2)
  }
}

cat("*** cargo is ok\n")
quit("no", status = 0)
