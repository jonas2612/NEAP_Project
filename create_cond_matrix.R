library(dplyr)
library(stringr)
library(jsonlite)

meta <- readxl::read_excel(
  path = "C:\\Users\\jonas\\OneDrive\\Desktop\\Semester_4_Bio\\Masterpraktikum\\metadata_all_samples_new.xlsx",
  sheet = 1,
  col_types = "text",
  .name_repair = "minimal"
  )|>
    mutate(
      location = str_trim(as.character(location)),
      medication = str_trim(as.character(medication)),
      KO = str_trim(as.character(KO)),
      dose = str_trim(as.character(dose)),
      duration = str_trim(as.character(duration)),
      CAD = str_trim(as.character(CAD)),
      CVD = str_trim(as.character(CVD)),
      combinedCondition = str_trim(as.character(combined_condition)),
      symptomaticAtherosclerosis = str_trim(as.character(symptomatic_atherosclerosis))
    )

# Reproduce Step 1 exactly, per-column
sanitize_group_levels <- function(raw_values) {
  raw_values <- str_trim(as.character(raw_values))
  levels_raw <- unique(raw_values)
  data.frame(
    groupraw  = levels_raw,
    groupcoef = make.names(levels_raw, unique = TRUE),
    stringsAsFactors = FALSE
  )
}

# Reproduce Step 2 exactly
safe_label_from_contrast <- function(test_coef, reference_coef) {
  contrast_expression <- paste0(test_coef, "_-_", reference_coef)
  x <- gsub("[^[:alnum:]._-]", "_", contrast_expression)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

build_label_lookup <- function(meta, groupcol) {
  meta |>
    group_by(accession) |>
    group_modify(~ {
      gk <- sanitize_group_levels(.x[[groupcol]])
      combos <- expand.grid(test = gk$groupraw, reference = gk$groupraw, stringsAsFactors = FALSE) |>
        filter(test != reference)
      combos |>
        left_join(gk, by = c("test" = "groupraw")) |>
        rename(test_coef = groupcoef) |>
        left_join(gk, by = c("reference" = "groupraw")) |>
        rename(reference_coef = groupcoef) |>
        mutate(
          groupcol = groupcol,
          safe_label = safe_label_from_contrast(test_coef, reference_coef)
        )
    }) |>
    ungroup()
}

label_lookup <- bind_rows(
  build_label_lookup(meta, "combinedcondition"),
  build_label_lookup(meta, "medication"),
  build_label_lookup(meta, "KO"),
  build_label_lookup(meta, "dose"),
  build_label_lookup(meta, "duration"),
  build_label_lookup(meta, "CAD"),
  build_label_lookup(meta, "CAD"),
  build_label_lookup(meta, "CVD"),
  build_label_lookup(meta, "symptomaticAtherosclerosis"),
  build_label_lookup(meta, "location"),
)

#####miRNA#########
json <- fromJSON("C:\\Users\\jonas\\OneDrive\\Desktop\\Semester_4_Bio\\Masterpraktikum\\relevant_views_miRNA.json", simplifyVector = FALSE)

accession_index_cutoff <- function(parts) {
  i = 1
  for (x in parts) {
    if (!startsWith(x, "GSE")) {
      break
    }
    i <- i + 1
  }
  i
}

extract_all_strings <- function(x) {
  nms <- names(x)
  if (is.null(nms)) return(character(0))
  
  sapply(nms, function(n) {
    parts <- strsplit(n, "_")[[1]]
    if (length(parts) <= 2) return("") # make it comb GSE agnostic
    message(accession_index_cutoff(parts))
    paste(parts[accession_index_cutoff(parts):(length(parts) - 1)], collapse = "_")
  }, USE.NAMES = FALSE)
}

json_names <- unique(str_trim(extract_all_strings(json)))

tmp <- names(json)[1]
parts <- strsplit(tmp, "_")[[1]]
out <- label_lookup[label_lookup$accession==paste(parts[1:accession_index_cutoff(parts)-1], collapse="_") & label_lookup$safe_label==paste(parts[accession_index_cutoff(parts):(length(parts) - 1)], collapse = "_"), ]# make it comb GSE agnostic

for (name in names(json)[2:length(names(json))]) {
  parts <- strsplit(name, "_")[[1]]
  print(paste(parts[1:accession_index_cutoff(parts)-1], collapse="_"))
  out <- rbind(out, label_lookup[label_lookup$accession==paste(parts[1:accession_index_cutoff(parts)-1], collapse="_") & label_lookup$safe_label==paste(parts[accession_index_cutoff(parts):(length(parts) - 1)], collapse = "_"), ])
}
out

write.csv(out, "C:\\Users\\jonas\\OneDrive\\Desktop\\Semester_4_Bio\\Masterpraktikum\\condition_matches_miRNA.csv")


######mRNA############
json <- fromJSON("C:\\Users\\jonas\\OneDrive\\Desktop\\Semester_4_Bio\\Masterpraktikum\\relevant_views_mRNA.json", simplifyVector = FALSE)

extract_all_strings <- function(x) {
  nms <- names(x)
  if (is.null(nms)) return(character(0))
  
  sapply(nms, function(n) {
    parts <- strsplit(n, "_")[[1]]
    if (length(parts) <= 2) return("")
    paste(parts[2:(length(parts) - 1)], collapse = "_")
  }, USE.NAMES = FALSE)
}

json_names <- unique(str_trim(extract_all_strings(json)))

tmp <- names(json)[1]
parts <- strsplit(tmp, "_")[[1]]
out <- label_lookup[label_lookup$accession==parts[1] & label_lookup$safe_label==paste(parts[2:(length(parts) - 1)], collapse = "_"), ]

for (name in names(json)[2:length(names(json))]) {
  parts <- strsplit(name, "_")[[1]]
  out <- rbind(out, label_lookup[label_lookup$accession==parts[1] & label_lookup$safe_label==paste(parts[2:(length(parts) - 1)], collapse = "_"), ])
}
out

write.csv(out, "C:\\Users\\jonas\\OneDrive\\Desktop\\Semester_4_Bio\\Masterpraktikum\\condition_matches_mRNA.csv")
