#' @keywords internal
.has_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

#' @keywords internal
.bioc_optional_pkgs <- function() {
  c(
    "AUCell", "GSVA", "slingshot", "tradeSeq", "DESeq2", "edgeR",
    "limma", "clusterProfiler", "ReactomePA", "monocle", "monocle3"
  )
}

#' @keywords internal
.default_install_hint <- function(pkg) {
  if (pkg %in% .bioc_optional_pkgs()) {
    sprintf("Install with BiocManager::install('%s').", pkg)
  } else {
    sprintf("Install with install.packages('%s').", pkg)
  }
}

#' @keywords internal
.monocle3_install_hint <- function() {
  paste(
    "Monocle3 is optional and only required for trajectory workflows.",
    "You can install it separately with BiocManager::install('monocle3').",
    "If you do not need trajectory analysis, you can ignore this dependency."
  )
}

#' @keywords internal
.assert_pkg <- function(pkg, feature = NULL, install_hint = NULL) {
  if (.has_pkg(pkg)) {
    return(invisible(TRUE))
  }

  if (is.null(install_hint)) {
    install_hint <- .default_install_hint(pkg)
  }

  msg <- if (is.null(feature)) {
    sprintf("Optional package '%s' is required.", pkg)
  } else {
    sprintf("Optional package '%s' is required for %s.", pkg, feature)
  }
  stop(paste(msg, install_hint), call. = FALSE)
}

#' @keywords internal
.assert_monocle3 <- function(feature = "trajectory analysis") {
  .assert_pkg(
    pkg = "monocle3",
    feature = feature,
    install_hint = .monocle3_install_hint()
  )
}
