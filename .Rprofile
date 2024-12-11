if (requireNamespace("logger", quietly = TRUE)) {
  library(logger)
} else {
  message("Package `logger` not available.")
}

## This makes sure that R loads the workflowr package
## automatically, everytime the project is loaded
# if (requireNamespace("workflowr", quietly = TRUE)) {
#   message("Loading .Rprofile for the current workflowr project")
#   library("workflowr")
# } else {
#   message("workflowr package not installed, please run install.packages(\"workflowr\") to use the workflowr functions")
# }

options(
  formatR.width = 79,
  formatR.args.newline = TRUE,
  formatR.brace.newline = FALSE,
  formatR.indent = 2
)


# DISPLAY envvar set to xpra session
if (stringr::str_detect(Sys.info()["nodename"], "rkn")) {
  Sys.setenv("DISPLAY" = ":103")
}

mkDir <- function(path) {
  if (!dir.exists(here(path))) {
    dir.create(here(path), recursive = TRUE)
  }
}


load_utils <- function() {
  util <- new.env()
  source(here("analysis/99_util.R"), local = util)
  print(ls.str(envir = util))
  return(util)
}


session_info <- function() {
  session <- sessionInfo()
  saveRDS(session, here(str_glue("logs/R/{prefix}_session.rds")))
  print(session)
}


save_versions <- function(...) {
  packages <- list(...)
  package_versions <- purrr::map(
    packages,
    (\(pkg) as.character(packageVersion(pkg)))
  ) |>
    purrr::set_names(packages)
  package_versions$R <- as.character(package_version(R.Version()))

  yaml_path <- here::here(stringr::str_glue("logs/R/{prefix}_versions.yaml"))
  yaml::write_yaml(package_versions, yaml_path)
  log_info("[VERSIONS] Saved package versions to {yaml_path}")
}


saveplot <- function(plotObj, path = NULL, h = 7, w = 7, name = NULL, ...) {
  objname <- ifelse(is.null(name), deparse(substitute(plotObj)), name)

  if (is.null(path)) {
    fileprefix <- here(str_glue("plots/{prefix}"))
    mkDir(fileprefix)
    filepath <- str_glue("{fileprefix}/{objname}.pdf")
  } else {
    filepath <- here(path)
  }

  if ("gg" %in% class(plotObj)) {
    tryCatch(
      {
        ggplot2::ggsave(
          filepath,
          width = w,
          height = h,
          plot = plotObj,
          ...
        )
        message(stringr::str_glue("Saved {objname} to {filepath}"))
        log_info("[PLOT]: Saved {objname} to {filepath}")
      },
      error = function(e) {
        print(e)
        log_warn("[PLOT]: {objname} might still be locked by git annex.")
      }
    )
  } else if ("Heatmap" %in% class(plotObj)) {
    tryCatch(
      {
        pdf(filepath,
          w = w, h = h
        )
        draw(plotObj)
        dev.off()
        message(stringr::str_glue("Saved {objname} to {filepath}"))
        log_info("[PLOT]: Saved {objname} to {filepath}")
      },
      error = function(e) {
        print(e)
        log_warn("[PLOT]: {objname} is probably still locked by git annex")
      }
    )
  } else {
    log_error("[PLOT] SOMETHING WENT WRONG?!?!")
  }
}




saveobj <- function(obj, path = NULL) {
  objname <- deparse(substitute(obj))

  if (is.null(path)) {
    fileprefix <- here(str_glue("results/R/{prefix}"))
    mkDir(fileprefix)
    filepath <- str_glue("{fileprefix}/{objname}.rds")
  } else {
    filepath <- here(path)
  }
  tryCatch(
    {
      saveRDS(obj, file = filepath)
      message(stringr::str_glue("Saved {objname} to {filepath}"))
      log_info("[OBJECT]: Saved {objname} to {filepath}")
    },
    error = function(e) {
      log_warn("[OBJECT]: {objname} is probably still locked by git annex")
    }
  )
}


objsize <- function(obj) {
  utils:::format.object_size(object.size(obj), "auto")
}


parseArguments <- function(arguments, description) {
  parser <- ArgumentParser(description = description)

  # Add Arugment groups
  map(toupper(names(arguments)), function(x) parser$add_argument_group(x))

  # Add arguments
  for (group in arguments) {
    map(group, function(x) do.call(parser$add_argument, x))
  }

  return(parser)
}


parseSnakefile <- function(rule) {
  code <- paste(
    "import snakemake",
    "import os",
    paste0("os.chdir('", here::here(), "')"),
    "workflow = snakemake.Workflow(snakefile = 'workflow/Snakefile', rerun_triggers = 'mtime')",
    "workflow.include('workflow/Snakefile')",
    paste0("rule_obj = workflow.get_rule('", rule, "')"),
    "input = {name: str(file) for name, file in rule_obj._input.items()}",
    "output = {name: str(file) for name, file in rule_obj._output.items()}",
    sep = "\n"
  )
  py <- reticulate::py_run_string(code)

  snakemake <- list(
    input = py$input,
    output = py$output
  )

  snakemake
}


colours <- list(
  studies = c(
    "Barzadd" = "#999999", "Chiang" = "#E69F00",
    "Dhiman" = "#56B4E9", "Hefzi" = "#009E73", "Kol" = "#F0E442",
    "Orellana" = "#0072B2", "Tzani" = "#D55E00", "van Wijk" = "#CC79A7"
  ),
  cellLines = c("CHO-S" = "red", "CHO-K1" = "black"),
  producers = c("TRUE" = "orange", "FALSE" = "white"),
  libraries = c("PolyA" = "#b9cd4f", "rRNA" = "#0191c1")
)
