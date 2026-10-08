library(multiMiR)
library(dplyr)
library(tidyr)
library(readr)

mirnas_mouse <- c(
  'mmu-miR-139-5p',
  'mmu-miR-3960',
  'mmu-miR-133a-3p',
  'mmu-miR-136-5p',
  'mmu-miR-21a-5p',
  'mmu-miR-146b-5p',
  'mmu-miR-23b-3p',
  'mmu-miR-15b-5p',
  'mmu-miR-7219-3p',
  'mmu-miR-101c',
  'mmu-miR-133b-3p',
  'mmu-miR-221-3p',
  'mmu-miR-143-5p',
  'mmu-miR-101a-3p',
  'mmu-let-7i-5p',
  'mmu-miR-714',
  'mmu-miR-224-5p',
  'mmu-miR-532-3p',
  'mmu-miR-222-3p',
  'mmu-miR-328-5p',
  'mmu-miR-26b-5p',
  'mmu-miR-146a-5p',
  'mmu-miR-26a-5p',
  'mmu-miR-301a-3p',
  'mmu-miR-342-5p',
  'mmu-miR-329-3p',
  'mmu-miR-34a-5p',
  'mmu-miR-342-3p',
  'mmu-miR-1a-3p',
  'mmu-miR-143-3p',
  'mmu-let-7e-3p',
  'mmu-miR-30b-5p',
  'mmu-miR-30d-5p',
  'mmu-miR-142a-3p'
)

res <- get_multimir(
  org="mmu", mirna=mirnas_mouse, table = "validated", summary=T
)
validated_tbl <- res@data %>% as_tibble()

res <- get_multimir(
  org="mmu", mirna=mirnas_mouse, table = "predicted", summary=T
)
predicted_tbl <- res@data %>% as_tibble()

val_pairs <- validated_tbl %>% transmute(
  miRNA=mature_mirna_id, gene=target_symbol, evidence="validated", source = database
)

pred_pairs <- predicted_tbl %>% transmute(
  miRNA=mature_mirna_id, gene=target_symbol, evidence="predicted", source = paste(sort(unique(database)), collapse = "; ")
)

combined <- bind_rows(val_pairs, pred_pairs) %>% filter(!is.na(gene), gene != "") %>% distinct()
combined

write_csv(combined, "C:\\Users\\jonas\\OneDrive\\Desktop\\Semester_4_Bio\\Masterpraktikum\\miRNA_targets_mouse.csv")



mirnas_human <- c(
  'hsa-miR-1260b',
  'hsa-miR-8069',
  'hsa-miR-494-3p',
  'hsa-miR-200b-3p',
  'hsa-miR-6085',
  'hsa-miR-20a-5p',
  'hsa-miR-100-5p',
  'hsa-miR-128-3p',
  'hsa-miR-17-3p',
  'hsa-let-7a-5p',
  'hsa-miR-5739',
  'hsa-miR-7977',
  'hsa-miR-30a-5p',
  'hsa-miR-5100',
  'hsa-miR-19b-3p',
  'hsa-let-7b-5p',
  'hsa-miR-4286',
  'hsa-miR-7975',
  'hsa-miR-92a-3p',
  'hsa-miR-6826-5p',
  'hsa-miR-1260a',
  'hsa-miR-6089',
  'hsa-miR-16-2-3p',
  'hsa-miR-4697-5p',
  'hsa-miR-1225-5p',
  'hsa-miR-4687-3p',
  'hsa-miR-6090',
  'hsa-miR-4516',
  'hsa-miR-6749-5p',
  'hsa-miR-1202',
  'hsa-miR-8072',
  'hsa-miR-6165',
  'hsa-miR-320b',
  'hsa-let-7f-5p',
  'hsa-miR-29b-3p'
)

res <- get_multimir(
  org="hsa", mirna=mirnas_human, table = "validated", summary=T
)
validated_tbl <- res@data %>% as_tibble()

res <- get_multimir(
  org="hsa", mirna=mirnas_human, table = "predicted", summary=T
)
predicted_tbl <- res@data %>% as_tibble()

val_pairs <- validated_tbl %>% transmute(
  miRNA=mature_mirna_id, gene=target_symbol, evidence="validated", source = database
)

pred_pairs <- predicted_tbl %>% transmute(
  miRNA=mature_mirna_id, gene=target_symbol, evidence="predicted", source = paste(sort(unique(database)), collapse = "; ")
)

combined <- bind_rows(val_pairs, pred_pairs) %>% filter(!is.na(gene), gene != "") %>% distinct()
combined

write_csv(combined, "C:\\Users\\jonas\\OneDrive\\Desktop\\Semester_4_Bio\\Masterpraktikum\\miRNA_targets_human.csv")


mirnas_human_short <- c(
  'hsa-miR-1260b',
  'hsa-miR-8069',
  'hsa-miR-494-3p',
  'hsa-miR-200b-3p',
  'hsa-miR-6085',
  'hsa-miR-20a-5p',
  'hsa-miR-100-5p',
  'hsa-miR-128-3p',
  'hsa-miR-17-3p',
  'hsa-let-7a-5p',
  'hsa-miR-5739',
  'hsa-miR-7977',
  'hsa-miR-30a-5p',
  'hsa-miR-5100',
  'hsa-miR-19b-3p',
  'hsa-let-7b-5p',
  'hsa-miR-4286',
  'hsa-miR-7975',
  'hsa-miR-92a-3p',
  'hsa-miR-6826-5p',
  'hsa-miR-1260a',
  'hsa-miR-6089',
  'hsa-miR-16-2-3p'
)

res <- get_multimir(
  org="hsa", mirna=mirnas_human_short, table = "validated", summary=T
)
validated_tbl <- res@data %>% as_tibble()

res <- get_multimir(
  org="hsa", mirna=mirnas_human_short, table = "predicted", summary=T
)
predicted_tbl <- res@data %>% as_tibble()

val_pairs <- validated_tbl %>% transmute(
  miRNA=mature_mirna_id, gene=target_symbol, evidence="validated", source = database
)

pred_pairs <- predicted_tbl %>% transmute(
  miRNA=mature_mirna_id, gene=target_symbol, evidence="predicted", source = paste(sort(unique(database)), collapse = "; ")
)

combined <- bind_rows(val_pairs, pred_pairs) %>% filter(!is.na(gene), gene != "") %>% distinct()
combined

write_csv(combined, "C:\\Users\\jonas\\OneDrive\\Desktop\\Semester_4_Bio\\Masterpraktikum\\miRNA_targets_human_short.csv")