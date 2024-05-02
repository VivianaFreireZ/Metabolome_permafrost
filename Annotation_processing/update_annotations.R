here::i_am("Annotation_processing/update_annotations.R")
library(here)
library(tidyverse)

dram_annotations_directory <- here("data", "DRAM_annotations_v4")
output_annotations_directory <- here("data", "DRAM_annotations_v4")

read_dram_annotations <- function(filename = "all_annotations.tsv") {
  annotations <- read_tsv(here(dram_annotations_directory, filename)) %>%
    rename(gene_id = 1)
}

dram_annotations <- read_dram_annotations(filename = "all_annotations_raw.tsv")

graftm_assignments <- tribble(
  ~graftm_call, ~gene_name,
  # dsrA
  "Root; dsrA; Acidobacteriota", "dsrA",
  "Root; dsrA; Acidobacteriota_2", "dsrA",
  "Root; dsrA; Acidobacteriota_3", "dsrA",
  "Root; dsrA; Bacteroidota", "dsrA",
  "Root; dsrA; Chloroflexota", "dsrA",
  "Root; dsrA; Desulfovibrionales", "dsrA",
  "Root; dsrA; Myxococcota", "dsrA",
  "Root; dsrA; Nitrospirota", "dsrA",
  "Root; dsrA; Syntrophia", "dsrA",
  "Root; dsrA; Syntrophia_2", "dsrA",
  "Root; dsrA; Syntrophobacteria", "dsrA",
  "Root; dsrA; Syntrophorhabdia", "dsrA",
  "Root; dsrA; Verrucomicrobiota", "dsrA",
  "Root; rdsrA; Alphaproteobacteria", "rdsrA",
  "Root; rdsrA; Gammaproteobacteria", "rdsrA",
  "Root; rdsrA; Myxococcota_2", "rdsrA",
  # dsrB
  "Root; dsrB", "dsrB",
  "Root; dsrB; Acidobacteriota", "dsrB",
  "Root; dsrB; Acidobacteriota_2", "dsrB",
  "Root; dsrB; Bacteroidota", "dsrB",
  "Root; dsrB; Chloroflexota", "dsrB",
  "Root; dsrB; Desulfovibrionales", "dsrB",
  "Root; dsrB; Myxococcota", "dsrB",
  "Root; dsrB; Nitrospirota", "dsrB",
  "Root; dsrB; Syntrophia", "dsrB",
  "Root; dsrB; Syntrophia_2", "dsrB",
  "Root; dsrB; Syntrophobacteria", "dsrB",
  "Root; dsrB; Syntrophorhabdia", "dsrB",
  "Root; dsrB; Verrucomicrobiota", "dsrB",
  "Root; rdsrB", "rdsrB",
  "Root; rdsrB; Alphaproteobacteria", "rdsrB",
  "Root; rdsrB; Gammaproteobacteria_2", "rdsrB"
)

new_ko_ids <- tribble(
  ~gene_name, ~new_ko_id,
  # dsrA
  "dsrA", "K111800",
  "rdsrA", "K111801",
  # dsrB
  "dsrB", "K111810",
  "rdsrB", "K111811"
)

graftm_inputs <- tribble(
  ~gene, ~ko_id,
  "dsrA", "K11180",
  "dsrB", "K11181"
  ) %>%
  mutate(
    data = map(gene, ~ read_tsv(here("Annotation_processing", ., "emerge_all_read_tax.tsv"), col_names = c("gene_id", "graftm_call")))
  ) %>%
  unnest(data) %>%
  left_join(graftm_assignments) %>%
  left_join(new_ko_ids)

# Annotations ko_id: change to KO from graftm_inputs new_ko_id
new_annotations <- dram_annotations %>%
  left_join(graftm_inputs %>% select(gene_id, new_ko_id)) %>%
  mutate(
    ko_id = map2_chr(ko_id, new_ko_id, ~ ifelse(is.na(.y), .x, .y))
  )

write_tsv(
  new_annotations %>% select(gene_id, ko_id),
  here(output_annotations_directory, "updated_ko_ids.tsv"),
  na = ""
  )

# Replace annotations ko_id column with new column
# DRAM_DIR=data/DRAM_annotations_v4
# awk 'BEGIN{OFS=FS="\t"}FNR==NR{a[NR]=$2;next}{$10=a[FNR]}1' $DRAM_DIR/updated_ko_ids.tsv $DRAM_DIR/all_annotations_raw.tsv > $DRAM_DIR/all_annotations.tsv
