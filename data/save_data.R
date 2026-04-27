###############################################################################
# Convert MATLAB replication data files to .rds format
# Source: ~/Downloads/DGU_Replication_Code_and_Data/Data/
###############################################################################

SRC_DIR <- "~/Downloads/DGU_Replication_Code_and_Data/Data"
OUT_DIR <- "data"

# Dataset definitions: name, source file, start_date filter, nFactors
datasets <- list(
  SPSectors = list(
    file       = "SPSectors.txt",
    start_date = NULL,
    nFactors   = 1
  ),
  Industry = list(
    file       = "Industry.txt",
    start_date = 19630728,
    nFactors   = 1
  ),
  International = list(
    file       = "international.txt",
    start_date = NULL,
    nFactors   = 1
  ),
  MKT_SMB_HML = list(
    file       = "F-F_Research_Data_Factors.txt",
    start_date = 19630728,
    nFactors   = 1
  ),
  FF_1factor = list(
    file       = "25_Portfolios_5x5_1Factor.txt",
    start_date = 19630728,
    nFactors   = 1
  ),
  FF_4factor = list(
    file       = "25_Portfolios_5x5_MOM.txt",
    start_date = 19630728,
    nFactors   = 4
  )
)

for (ds_name in names(datasets)) {
  cfg <- datasets[[ds_name]]
  filepath <- file.path(SRC_DIR, cfg$file)

  # Read with read.table first; fallback to read.delim for files with empty fields
  raw <- tryCatch(
    read.table(filepath, header = FALSE, skip = 1, stringsAsFactors = FALSE),
    error = function(e) {
      read.delim(filepath, header = FALSE, skip = 1, sep = "\t",
                 fill = TRUE, stringsAsFactors = FALSE)
    }
  )

  # Filter to start date
  if (!is.null(cfg$start_date)) {
    row_start <- which(raw[, 1] == cfg$start_date)
    if (length(row_start) == 0) stop("Start date not found for ", ds_name)
    raw <- raw[row_start:nrow(raw), ]
  }

  # Remove incomplete rows
  raw <- raw[complete.cases(raw), ]

  # Read header for column names
  header_line <- readLines(filepath, n = 1)
  header_line <- sub("^%", "", header_line)
  col_names <- trimws(strsplit(header_line, "\t")[[1]])
  # Handle whitespace-separated headers
  if (length(col_names) < ncol(raw)) {
    col_names <- trimws(strsplit(header_line, "\\s+")[[1]])
  }
  if (length(col_names) == ncol(raw)) {
    colnames(raw) <- col_names
  } else {
    colnames(raw) <- paste0("V", seq_len(ncol(raw)))
    colnames(raw)[1] <- "Date"
    colnames(raw)[2] <- "RF"
  }

  # Store metadata as attributes
  dat <- list(
    data     = raw,
    nFactors = cfg$nFactors,
    source   = cfg$file
  )

  out_path <- file.path(OUT_DIR, paste0(ds_name, ".rds"))
  saveRDS(dat, out_path)
  cat(sprintf("Saved %s: %d rows x %d cols -> %s\n",
              ds_name, nrow(raw), ncol(raw), out_path))
}
