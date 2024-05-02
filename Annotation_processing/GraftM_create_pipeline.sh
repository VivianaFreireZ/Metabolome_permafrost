#!/bin/bash
# environment: mmseqs / graftm

OUTPUT_DIR=results/
SEEDS=inputs/seeds.faa
SEEDS_TAX=inputs/seeds_tax.tsv
KO=K00001
GENE=genE
mkdir -p $OUTPUT_DIR

#############
### Setup ###
#############
UNIREF_DB=uniref/release_2022_01/uniref90
UNIREF_FASTA=uniref/release_2022_01/uniref90.fasta
EMERGE_GENES_FASTA=bins/combined_all_genes.faa
EMERGE_GENES_IDS=data/gene_ids_for_GraftM.tsv

##################################
### Search Uniref90 for $SEEDS ###
##################################
mmseqs easy-search \
    --threads 10 \
    $SEEDS \
    $UNIREF_DB \
    $OUTPUT_DIR/uniref90_search.m8 \
    $OUTPUT_DIR/tmp \
    &> $OUTPUT_DIR/mmseqs.log

# Extract hit sequences (300 sequences per seed sequence)
# Switch to GraftM env
awk -F "\t" '{print $2}' $OUTPUT_DIR/uniref90_search.m8 | sort | uniq > $OUTPUT_DIR/uniref90_search_ids.txt

mfqe \
    --sequence-name-lists $OUTPUT_DIR/uniref90_search_ids.txt \
    --input-fasta $UNIREF_FASTA \
    --output-fasta-files $OUTPUT_DIR/uniref90_search.faa \
    --output-uncompressed \
    &> $OUTPUT_DIR/mfqe.log

cat $SEEDS $OUTPUT_DIR/uniref90_search.faa > $OUTPUT_DIR/uniref90_search_and_seeds.faa

############################
### Create taxonomy file ###
############################
sed "s/$/\t${GENE}/" $OUTPUT_DIR/uniref90_search_ids.txt > $OUTPUT_DIR/uniref90_search_tax.txt

cat $SEEDS_TAX $OUTPUT_DIR/uniref90_search_tax.txt > $OUTPUT_DIR/uniref90_search_and_seeds_tax.txt

############################
### Add EMERGE sequences ###
############################
# ko_id, gene id, genome id (exported from Metabolic-analysis/exports.Rmd)
grep $KO $EMERGE_GENES_IDS > $OUTPUT_DIR/emerge_info.txt
awk -F "\t" '{print $2}' $OUTPUT_DIR/emerge_info.txt > $OUTPUT_DIR/emerge_ids.txt
mfqe \
    --sequence-name-lists $OUTPUT_DIR/emerge_ids.txt \
    --input-fasta $EMERGE_GENES_FASTA \
    --output-fasta-files $OUTPUT_DIR/emerge.faa \
    --output-uncompressed \
    &> $OUTPUT_DIR/mfqe_emerge.log

awk 'BEGIN{FS=OFS="\t"}{print $2, "emerge;" $3}' $OUTPUT_DIR/emerge_info.txt > $OUTPUT_DIR/emerge_tax.txt

#############################
### Create GraftM package ###
#############################
# Combine all sequences and taxonomy
cat $OUTPUT_DIR/uniref90_search_and_seeds.faa $OUTPUT_DIR/emerge.faa > $OUTPUT_DIR/combined.faa
cat $OUTPUT_DIR/uniref90_search_and_seeds_tax.txt $OUTPUT_DIR/emerge_tax.txt > $OUTPUT_DIR/combined_tax.txt

graftM create \
    --sequences $OUTPUT_DIR/combined.faa \
    --taxonomy $OUTPUT_DIR/combined_tax.txt \
    --output $OUTPUT_DIR/${GENE}_draft.gpkg \
    &> $OUTPUT_DIR/GraftM_draft.log

#-- If automated rooting fails
# Manual rerooting with arb at base of separate gene type
# Load alignment from $OUTPUT_DIR/${GENE}_draft.gpkg/${GENE}_draft.gpkg.refpkg/combined_deduplicated_aligned.fasta
# Load tree from $OUTPUT_DIR/${GENE}_draft.gpkg/${GENE}_draft.gpkg.refpkg/*.tree
# Export rerooted tree to $OUTPUT_DIR/combined_rerooted.tree
graftM create \
    --taxtastic_taxonomy graftm_create_taxonomy.combined.csv \
    --taxtastic_seqinfo graftm_create_seqinfo.combined.csv \
    --alignment graftm_create_alignment.combined.faa  \
    --rerooted_tree $OUTPUT_DIR/combined_rerooted.tree \
    --sequences $OUTPUT_DIR/combined.faa \
    --output $OUTPUT_DIR/${GENE}_draft.gpkg \
    &> $OUTPUT_DIR/GraftM_reroot.log
#-- End

##################################
### Extract sequences with Arb ###
##################################
# Load alignment from $OUTPUT_DIR/${GENE}_draft.gpkg/${GENE}_draft.gpkg.refpkg/combined_deduplicated_aligned.fasta
# Load tree from $OUTPUT_DIR/combined_rerooted.tree or  $OUTPUT_DIR/${GENE}_draft.gpkg/${GENE}_draft.gpkg.refpkg/tree*.tre
# Import calc-sheet from $OUTPUT_DIR/combined_tax.txt
## col 2 to field: taxonomy
## where name matches col 1
# Add taxonomy (100 width) to NDS in Tree (remove name, acc)

# Colour important sequences:
## normal seeds: 2
## reverse seeds: 1
## reverse UniRef: 3
## Emerge sequences: 4
# Group sequences into clades by taxonomy (hopefully defined by seed sequences)
# Compare tree to published trees

# Export tree in newick format with group names
## Filename: $OUTPUT_DIR/manual_tax.tree

#############################
### Remake GraftM package ###
#############################
# Create package based on annotated tree
graftM create \
    --sequences $OUTPUT_DIR/combined.faa \
    --alignment $OUTPUT_DIR/${GENE}_draft.gpkg/${GENE}_draft.gpkg.refpkg/combined_deduplicated_aligned.fasta \
    --rerooted_annotated_tree $OUTPUT_DIR/manual_tax.tree \
    --output $OUTPUT_DIR/$GENE.gpkg \
    &> $OUTPUT_DIR/GraftM.log

##########################
### Test final package ###
##########################
# Run against EMERGE ${GENE} KO sequences
graftM graft \
    --forward $OUTPUT_DIR/emerge.faa \
    --graftm_package $OUTPUT_DIR/$GENE.gpkg \
    --output_directory $OUTPUT_DIR/GraftM_output \
    &> $OUTPUT_DIR/GraftM_graft.log

# Result in $OUTPUT_DIR/GraftM_output/emerge/emerge_read_tax.tsv
