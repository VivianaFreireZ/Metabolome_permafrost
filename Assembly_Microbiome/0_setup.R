# Install packages
here::i_am("setup.R")
library(here)
library(tidyverse)
library(ape)
library(phytools)
library(picante)

# Set V1 directory
v1_mags_directory <- here("data", "Emerge_MAGs_v1")

read_genome_info <- function() {
    d <- bind_rows(
        read_tsv(here(v1_mags_directory, "Field_bins", "Woodcroft_2018", "tax_complete_contam.txt")) %>% mutate(genome_set = "Woodcroft_2018"),
        read_tsv(here(v1_mags_directory, "Field_bins", "Cronin_2021", "tax_complete_contam.txt")) %>% mutate(genome_set = "Cronin_2021"),
        read_tsv(here(v1_mags_directory, "SIP_bins", "tax_complete_contam.txt")) %>% mutate(genome_set = "SIP_bins"),
        read_tsv(here(v1_mags_directory, "JGI_bins", "tax_complete_contam.txt")) %>% mutate(genome_set = "JGI_bins"),
        )
    return(d)
}

# Load the MAG abundances
read_trimmed_mean <- function(trimmed_mean_file = "MAG_otu_Trimmed_Mean.txt",
                              SampleID_name = "hiseq_sample_id") {
    d <- read_csv(here("data", trimmed_mean_file)) %>%
        rename(genome = Genome) %>%
        pivot_longer(cols = !genome, names_to = SampleID_name, values_to = "coverage") %>%
        filter(coverage > 0) %>%
        group_by(!!as.name(SampleID_name)) %>%
        mutate(relabund_of_recovered = coverage / sum(coverage))
    
    rel_abund <- d %>%
        select(genome, !!as.name(SampleID_name), relabund_of_recovered) %>%
        pivot_wider(names_from = !!as.name(SampleID_name), values_from = relabund_of_recovered,
                    values_fill = 0)
    coverage <- d %>%
        select(genome, !!as.name(SampleID_name), coverage) %>%
        pivot_wider(names_from = !!as.name(SampleID_name), values_from = coverage,
                    values_fill = 0)
    
    read_trimmed_mean <- list(trimmed_mean = d, 
                              rel_abund = rel_abund, coverage = coverage)
    
    return(read_trimmed_mean)
}

# Use temporal_sample_id to match
read_sample_metadata <- function() {
    d <- read_csv(here("data", "coring_geochem_sequenced_samples.csv")) %>%
        rename(hash = 1) %>%
        mutate(
            sample = SampleID__,
            temporal_sample_id = metaG_JGI_NovaSeq__,
            temporal_sample_id = gsub(" ", "_", temporal_sample_id)
        )
    # Read in Vivi's metadata
    v <- read_csv(here("data", "metadata_matched_MAGs_data_11-5-21.csv"))

    # Merge the Microbiome Metadata sheet with Vivi's so that only her samples are 
    # retained
    d <- left_join(v, d, by = c("SampleID_matching_MAGs" = "SampleID__")) %>%
      mutate(SampleID_matching_MAGs = metaG_ACE__) %>%
      # fix the two missing samples manually
      mutate(SampleID_matching_MAGs = ifelse(SampleID == "July_S_3_M_12_so_R1", "20120700_S3M", SampleID_matching_MAGs),
             SampleID_matching_MAGs = ifelse(SampleID == "Aug_P_3_S_12_so_R1", "20120800_P3S", SampleID_matching_MAGs))

    return(d)
}

read_genome_tree <- function(otu_table, visualize = FALSE) {
    require(ape)
    require(phytools)
    
    arch.tree.raw <- phytools::read.newick(here("data", "gtdbtk.ar122.decorated.tree"))
    bac.tree.raw <- phytools::read.newick(here("data", "gtdbtk.bac120.decorated.tree"))
    
    # Artificially merge the two trees based on the rooting proposed in paper out of 
    # Rob Knight's lab: Zhu et al. 2019 https://www.nature.com/articles/s41467-019-13443-4: 
    # "Phylogenomics of 10,575 genomes reveals evolutionary proximity between domains Bacteria and Archaea"
    bac.arch.tree.raw <- ape::bind.tree(bac.tree.raw, arch.tree.raw)
    
    # check tree merge
    if(visualize) {
        writeLines("Plotting initial tree")
        plot(bac.arch.tree.raw)
    }
    
    # Fix the names by removing the extra quotes
    bac.arch.tree.raw$tip.label <- gsub("'", "", bac.arch.tree.raw$tip.label)
    
    # Drop tips that aren't from EMERGE MAGs 
    #(since we used gtdb-tk to make the tree these should be tips that are in gtdb, but not our MAGs)
    gtdb_tips <- bac.arch.tree.raw$tip.label[!(bac.arch.tree.raw$tip.label %in% otu_table$genome)]
    
    bac.arch.tree.nogtdb <- ape::drop.tip(bac.arch.tree.raw, tip = gtdb_tips)
    bac.tree <- ape::drop.tip(bac.tree.raw, tip = gtdb_tips)
    arch.tree <- ape::drop.tip(arch.tree.raw, tip = gtdb_tips)
    
    if(visualize) {
        writeLines("Plotting filtered tree")
        # plot resulting tree
        plot(bac.arch.tree.nogtdb)
    }
    
    # Root the tree
    ba.tree.rt <- phytools::midpoint.root(bac.arch.tree.nogtdb)
    
    # check tree
    if(visualize) {
        writeLines("Plotting filtered tree")
        # plot resulting tree
        plot(ba.tree.rt)
    }
    return(list(bacteria_tree = bac.tree, 
                archaea_tree = arch.tree,
                bac_arch_tree = ba.tree.rt))
}



# Check for matching MAGs between datasets and get taxonomy
get_taxonomy <- function(otu_table, genome_info) {
    # Find genomes missing from OTU table (if any)
    missing_genomes_data <- otu_table$genome[!(otu_table$genome %in% 
                                               genome_info$genome)]
    writeLines(paste("Genomes missing from the otu_table that are present in",
                     "genome_info:", paste(missing_genomes_data, collapse = ", ")))
 
    
    # Prepare taxonomy table
    taxonomy_table <- genome_info %>%
        filter(genome %in% otu_table$genome) %>% # select only genomes in current data
        select(genome, taxonomy) %>%
        separate(taxonomy, 
                 into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
                 sep = ";") %>%
        select(genome, Domain:Species)
    
    return(taxonomy_table)
}

clean_sample_metadata <- function(sample_metadata) {
    sample_metadata_mod <- sample_metadata %>%
        # Modify depth code so that it is consistent across years
        # S = top 5 cm; M = 5-10cm, D = 10-15 CM, 15-20 = X0, all other
        # depths increment by 10cm and are X1-X6 with X6 being the deepest
        mutate(DepthCodeMod = ifelse(DepthAvg__ >= 1 & DepthAvg__ < 5, "S",
                              ifelse(DepthAvg__ >= 5 & DepthAvg__ < 10, "M", 
                              ifelse(DepthAvg__ >= 10 & DepthAvg__ < 15, "D", 
                              ifelse(DepthAvg__ >= 15 & DepthAvg__ <= 20, "X0",
                              ifelse(DepthAvg__ > 20 & DepthAvg__ < 30, "X1",
                              ifelse(DepthAvg__ >= 30 & DepthAvg__ < 40, "X2",
                              ifelse(DepthAvg__ >= 40 & DepthAvg__ < 50, "X3",
                              ifelse(DepthAvg__ >= 50 & DepthAvg__ < 60, "X4",
                              ifelse(DepthAvg__ >= 60 & DepthAvg__ < 70, "X5",
                              ifelse(DepthAvg__ >= 60 & DepthAvg__ < 80, "X6",
                              NA))))))))))) %>%
        mutate(DepthCodeMod = factor(DepthCodeMod, 
                                     levels = c("S", "M", "D","X0","X1","X2","X3","X4","X5","X6"))) %>%
        # ALD and temp fill in NAs for samples that have values in 
        # other columns. This fixes a quirk with how the data was entered
        # in the original field data sheets in one year that meant only
        # the first sample for a core had a recorded ALD and Air temp.
        group_by(Year__, Month__, Habitat__, Site__, Core__) %>%
        fill(ALD.cm__) %>%
        fill(`T, air (C)`) %>%
        ungroup() %>%
        # Try to solve the WTD/ALD detection issue (namely that when WTD
        # is below detection, we want to record this b/c it's particularly
        # relevant biologically but that the limit of detection changes
        # from year to year as the instruments changed)
        # Set samples "below detection" to -101 (lower than all observed values)
        # for WTD and to 101 for ALD (lower than all observed values and lower 
        # than the longest probe used); All these saved in new variables (e.g. 
        # "ALD" rather than "ALD.cm__") so they don't change any code that is 
        # dependent on the old versions.
        mutate(WTD = ifelse(grepl("Detect", WTD.cm_neg.is.below.sfc__), "-101", 
                            WTD.cm_neg.is.below.sfc__),
               WTD = ifelse(grepl("notes|Difficult|Pending|pending", WTD.cm_neg.is.below.sfc__), 
                            NA, WTD),
               WTD = as.numeric(WTD)) %>% 
        mutate(ALD = ifelse(grepl("Detect", ALD.cm__), "101", 
                            ALD.cm__),
               ALD = as.numeric(ALD)) %>% 
        # Set the measurements we're unsure of being less than or equal to 0 to NA
        mutate(`T, soil (C)` = ifelse(grepl("<=0", `T, soil (C)`), NA, `T, soil (C)`),
               `T, soil (C)` = as.numeric(`T, soil (C)`))
    
    return(sample_metadata_mod)
}


# Check for matching names between datasets
check_sample_metadata <- function(otu_table, sample_metadata, 
                                  sample_id_column = "SampleID_matching_MAGs",
                                  taxonomy = NULL,
                                  tree = NULL) {
    # Checking OTU table and sample metadata
    # Find the missing samples for otu_table and sample metadata
    missing_genomic_data <- na.omit(sample_metadata[[sample_id_column]][!(sample_metadata[[sample_id_column]] %in% names(otu_table[-1]))]) # Samples missing from OTU table
    writeLines(paste(length(missing_genomic_data)," samples are missing from the OTU table that are present in the",
                     "metadata:", paste(missing_genomic_data, collapse = ", ")))
    
    missing_env_data <- names(otu_table[-1])[!(names(otu_table[-1]) %in% sample_metadata[[sample_id_column]] )] # Samples missing from env data
    writeLines(paste(length(missing_env_data)," samples are  missing from the metadata that are present in the",
                     "OTU table:", paste(missing_env_data, collapse = ", ")))
    writeLines(paste0("Now filtering out missing sample(s)..."))
    
    # Filter out the samples that don't occur in both
    otu_table <- otu_table %>% select(genome, !all_of(missing_env_data))
    
    sample_metadata <- sample_metadata %>% 
        filter((!!as.name(sample_id_column) %in% names(otu_table[-1])))
    
    # After filtering samples, remove any columns which are entirely NAs (or 0s for MAGs)
    sample_metadata <- sample_metadata[,colSums(is.na(sample_metadata))<nrow(sample_metadata)]
    
    # Identify env_data columns that are all NA
    all_na_cols <- names(sample_metadata[,colSums(is.na(sample_metadata))==nrow(sample_metadata)])
    # filter out metadata variables that are all NAs
    sample_metadata <- sample_metadata[,colSums(is.na(sample_metadata))<nrow(sample_metadata)]

    
    # Reorder otu table columns to match order in sample metadata
    otu_table <- otu_table[,c("genome", sample_metadata[[sample_id_column]])]

    writeLines(paste(nrow(sample_metadata), "samples in sample_metadata"))
    writeLines(paste(ncol(otu_table[-1]), "samples in otu_table"))
    

    # Checking taxonomy, tree, and otu table (if taxonomy and tree are present)
    # and reordering to match tree order
    if(!is.null(tree) & !is.null(taxonomy)) {
        # reorder otu table to match order of tips in tree
        match.phylo.otu <- picante::match.phylo.data(tree, 
                                                     otu_table %>%
                                                         column_to_rownames(var = "genome"))
        tree <- match.phylo.otu$phy
        otu_table <- match.phylo.otu$data %>%
            rownames_to_column(var = "genome") %>%
            mutate(rowname = genome) %>%
            column_to_rownames()
        
        # Check that otu table genomes are all in the tree and taxonomy file
        genomes_missing_from_taxonomy <- otu_table$genome[!otu_table$genome %in% taxonomy$genome]
        writeLines(paste("Genomes missing from the taxonomy that are present in the",
                         "OTU table:", paste(genomes_missing_from_taxonomy, collapse = ", ")))
        
        genomes_missing_from_tree <- otu_table$genome[!otu_table$genome %in% tree$tip.label]
        writeLines(paste("Genomes missing from the tree that are present in the",
                         "OTU table:", paste(genomes_missing_from_tree, collapse = ", ")))
        
        writeLines(paste0(nrow(otu_table), " taxa in otu table"))
        writeLines(paste0(length(tree$tip.label), " tips in tree"))
        writeLines(paste0(nrow(taxonomy), " taxa in taxonomy"))
                                                   
        # reorder taxonomy to match otu_table and tree tip order
        taxonomy <- taxonomy %>% column_to_rownames("genome")
        taxonomy <- taxonomy[rownames(match.phylo.otu$data),] %>%
            rownames_to_column(var = "genome") %>%
            mutate(rowname = genome) %>%
            column_to_rownames()

    }

    
    # Create list of outputs to return:
    return_list <- list(sample_metadata = sample_metadata, otu_table = otu_table)
    
    if(!is.null(taxonomy)) {
        return_list$taxonomy <- taxonomy
    }
    
    if(!is.null(tree)) {
        return_list$tree <- tree
    }
    
    return(return_list)
}

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################

genome_info = read_genome_info()
trimmed_mean = read_trimmed_mean(trimmed_mean_file = "Trimmed_Mean_hiseq.txt")
sample_metadata = read_sample_metadata()

sample_metadata <- clean_sample_metadata(sample_metadata = sample_metadata)
trees <- read_genome_tree(otu_table = trimmed_mean$rel_abund,
                         visualize = FALSE)
taxonomy <- get_taxonomy(otu_table = trimmed_mean$rel_abund, 
                          genome_info = genome_info)



# Check sample metadata and combine into one output
input_ra <- check_sample_metadata(otu_table = trimmed_mean$rel_abund,
                      sample_metadata = sample_metadata,
                      taxonomy = taxonomy,
                      tree = trees$bac_arch_tree)

input_counts <- check_sample_metadata(otu_table = trimmed_mean$coverage,
                                  sample_metadata = sample_metadata,
                                  taxonomy = taxonomy,
                                  tree = trees$bac_arch_tree)
