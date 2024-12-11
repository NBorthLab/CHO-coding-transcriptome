suppressPackageStartupMessages({
    library(tidyverse)
    library(here)
    library(rio)
    library(argparse)
    # library(SRAdb)
    library(logger)
})


#' HIERARCHY OF NCBI SRA DATA
#'
#'        PROJECT
#'           |
#'        SAMPLE
#'           |
#'       EXPERIMENT
#'           |
#'          RUN


# ================================================
#   Parse command line
# ================================================

parser <- ArgumentParser()

parser$add_argument("--output",
    help = "Output in csv format",
    required = TRUE
)
# parser$add_argument(
#     "--db-file",
#     help = "Location of SRAdb sqlite file",
#     required = TRUE
# )
parser$add_argument("projects",
    nargs = "+",
    help = "List of Project numbers."
)

if (interactive()) {
    projects <- dir(here("resources/sra_runinfo"),
                    pattern = "_runinfo.csv",
                    full.names = TRUE)
    args <- parser$parse_args(c(
        "--output", "resources/sra_runinfo/all_runs.csv",
        projects
    ))
} else {
    args <- parser$parse_args()
}

log_info("Running with the following arguments")
print(args)

# ================================================
#   SRA databse
# ================================================


# NOTE: BROKEN AS FUCK
# Download runinfo with bash utilities
# See: https://github.com/seandavi/SRAdb/issues/37
# if (!file.exists(args$db_file)) {
#     getSRAdbFile(
#         destdir = here(dirname(args$db_file)),
#         destfile = basename(args$db_file)
#     )
# } else {
#     print("Using existing database file!")
# }


# ================================================
#   TSV input & Runs
# ================================================


# Import all CSV runinfos
runinfos <- args$projects %>%
    purrr::map_dfr(function(f) {
        rio::import(here(f)) %>%
            as_tibble() %>%
            # Library Name column should be a chracter
            mutate(LibraryName = as.character(LibraryName),
                   Sex = as.character(Sex)) %>%
            # Use only the transcriptomic RNA-seq runs
            dplyr::filter(LibraryStrategy == "RNA-Seq")
    })


# ================================================
#   Metadata
# ================================================

# trimmed_datasets <- c("")
# truseq_datasets <- c("SRP069883", "ERP122753", "SRP066848", "SRP111368")
#
runinfos <-
    runinfos %>%
    mutate(
        Study = case_when(
            SRAStudy == "SRP246348" ~ "Kol",
            SRAStudy == "SRP066848" ~ "van Wijk",
            SRAStudy == "SRP111368" ~ "Orellana",
            SRAStudy == "SRP159459" ~ "Chiang",
            SRAStudy == "SRP069883" ~ "Hefzi",
            SRAStudy == "SRP234382" ~ "Tzani",
            SRAStudy == "SRP324587" ~ "Barzadd",
            SRAStudy == "ERP122753" ~ "Dhiman"
        ),
        cell_lines = case_when(
            Study %in% c("Hefzi", "Kol") ~ "CHO-S",
            TRUE                         ~ "CHO-K1"
        ),
        producer = case_when(
            Study %in% c("van Wijk", "Chiang", "Hefzi") ~ FALSE,
            TRUE                                        ~ TRUE
        ),
        library_strat = case_when(
            Study %in% c("Dhiman", "Tzani") ~ "rRNA",
            TRUE                            ~ "PolyA"
        )
    )


# ================================================
#   Export
# ================================================

runinfos %>% rio::export(here(args$output))

log_info("Exported to {args$output}")
cat("\n")


# ================================================
#   Expermiments
# ================================================

# Get Runs of individual Experiments
# runinfos %>%
#     group_split(Experiment) %>%
#     map(function(x) {
#         if (nrow(x) > 1) {
#             x$Run
#         }
#     })
