#!/usr/bin/env Rscript --vanilla
# PRE RELEASE VERSION 
# EARLY ACCESS VERSION with NEW FEATURES  such as 
# BLAST analysis of the core genes for functional enrichment
# MetagenomeIntegrated has been relabelled as MetaPanR
# #####################################################################
# ################### Usage of MetaPanRan #################
# ####################################################################
# (base) git:(main) ✗ ./MetaPanRan.R -h
# usage: metagenomics_Integrated_version5.R [--] [--help] [--opts OPTS]
# [--annoatate_genome ANNOATATE_GENOME] [--roary_run ROARY_RUN]
# [--pEppan_run PEPPAN_RUN] [--pIrate_run PIRATE_RUN]
# [--pAnroo_run PANROO_RUN] [--gff_dir GFF_DIR] [--metadata_path
#  METADATA_PATH] [--cluster_meta_path CLUSTER_META_PATH]
# [--fasta_dir FASTA_DIR] [--align ALIGN] [--venn_viz VENN_VIZ]
# [--blast BLAST]
# 
# An integrated pipeline for the prediction of the bacterial genomes and
# analysis and visualization of the pangenomes from ROARY, PIRATE,
# PEPPAN, PANROO. This allows you to run the three accurate pangenome
# predictions from ROARY, PEPPAN, PIRATE, PANROO and then gives the
# pangenomes. It also gives a comparative assessment of the pangenomes
# that is conserved across all through a venn visualization. It also
# allows you to do the BLAST searches and also allows you to link the
# expression analysis of the core genes
# 
# flags:
#   -h, --help               show this help message and exit
# 
# optional arguments:
#   -x, --opts               RDS file containing argument values
# -a,                      --annoatate_genome. If you have just 
#                             assembled your genome your can call 
#                             this option to annotate the genome
# -r, --roary_run          Run roary as a part of the analysis
# -p, --pEppan_run         Will run the PEPPAN analysis and will
#                          convert to the roary for pangenome analysis
# --pIrate_run             Will run the PIRATE analysis and will
#                          Convert from PIRATE output to roary for
#                          pangenome analysis
# --pAnroo_run             Will run the PIRATE analysis and will
#                          Convert from PIRATE output to roary for
#                          pangenome analysis
# -g, --gff_dir            PATH to the gff directory to look for the
#                          presenece and absence gene
# -m, --metadata_path      add metadata for the pangenome
# -c, --cluster_meta_path  path for the meta data for the clusters
# -f, --fasta_dir          Provide the path to the folder containing
#                          the fasta file. You can provide the path to
#                          the directory
# --align                  align the selected genes using the IQTREE
#                          and the fastTREE using the MFP+MERGE model
#                          by finding the best partition scheme
#                          followed by tree inference and bootstrap and
#                          will also do a outlier detection in the
#                          aligned sequences and will plot a phylogeny
# -v, --venn_viz           create a venn visualiztion of the common
#                          pangenome across all the roary, pirate,
#                          peppan, panroo, For this you have to provide
#                          the presence and absence file from each of
#                          the run after running the entire workflow or
#                          you can provide directly the file.
# -b, --blast              this will blast all the genes selected from
#                          the pangenome for the similarity and also
#                          will prepare the blast files for the
#                          analysis
# ####################################################################
# ################### Arguments ######################################
# ####################################################################

suppressPackageStartupMessages(library(argparser, pos = "package:base"))
suppressPackageStartupMessages(library(methods, pos = "package:base"))
suppressPackageStartupMessages(library(pagoo, pos = "package:base"))
suppressPackageStartupMessages(library(Biostrings, pos = "package:base"))
suppressPackageStartupMessages(library(phangorn, pos = "package:base"))
suppressPackageStartupMessages(library(ggmsa, pos = "package:base"))
suppressPackageStartupMessages(library(ggplot2, pos = "package:base"))
suppressPackageStartupMessages(library(ape, pos = "package:base"))
suppressPackageStartupMessages(library(odseq, pos = "package:base"))
suppressPackageStartupMessages(library(PanVizGenerator, pos = "package:base"))
suppressPackageStartupMessages(library(ggmsa, pos = "package:base"))
suppressPackageStartupMessages(library(venn, pos = "package:base"))
suppressPackageStartupMessages(library(msa, pos = "package:base"))
suppressPackageStartupMessages(library(DECIPHER, pos = "package:base"))
suppressPackageStartupMessages(library(ggtree, pos = "package:base"))
suppressPackageStartupMessages(library(phytools, pos = "package:base"))
suppressPackageStartupMessages(library(reticulate, pos = "package:base"))

p <- arg_parser(
 "  An integrated pipeline for the prediction of the bacterial genomes and 
    analysis and visualization of the pangenomes from ROARY, PIRATE, PEPPAN, 
    PANROO. This allows you to run the three accurate pangenome predictions 
    from ROARY, PEPPAN, PIRATE, PANROO and then gives the pangenomes. It also 
    gives a comparative assessment of the pangenomes that is conserved across 
    all through a venn visualization. It also allows you to do the BLAST searches
    and also allows you to link the expression analysis of the core genes"
)
p <- add_argument(p, "--annoatate_genome", 
                  help = "If you have just assembled your genome
                  your can call this option to annotate the genome")
p <- add_argument(p, "--roary_run",
                  help = "Run roary as a part of the analysis")
p <- add_argument(p, "--pEppan_run",
                  help = "Will run the PEPPAN analysis and
                       will convert to the roary for pangenome 
                       analysis")
p <- add_argument(p, "--pIrate_run",
                  help = "Will run the PIRATE analysis and 
                  will Convert from PIRATE output to roary
                  for pangenome analysis")
p <- add_argument(p, "--pAnroo_run",
                  help = "Will run the PIRATE analysis and 
                  will Convert from PIRATE output to roary
                  for pangenome analysis")
p <- add_argument(p, "--gff_dir",
                  help = "PATH to the gff directory to look
                       for the presenece and absence gene")
p <- add_argument(p, "--metadata_path",
                  help = "add metadata for the pangenome")
p <- add_argument(p, "--cluster_meta_path",
                  help = "path for the meta data for the
                  clusters")
p <- add_argument(p,"--fasta_dir",
                  help = "Provide the path to the folder
                      containing the fasta file. You can
                      provide the path to the directory")
p <- add_argument(p, "--align",
                  help = "align the selected genes using
                  the IQTREE and the fastTREE using the 
                  MFP+MERGE model by finding the best 
                  partition scheme followed by tree 
                  inference and bootstrap and will also
                  do a outlier detection in the aligned
                  sequences and will plot a phylogeny")
p <- add_argument(p, "--venn_viz", 
                  help = "create a venn visualiztion of 
                  the common pangenome across all the 
                  roary, pirate, peppan, panroo, For this
                  you have to provide the presence and
                  absence file from each of the run
                  after running the entire workflow
                  or you can provide directly the file.")
p <- add_argument(p, "--blast", 
                  help = "this will blast all the genes 
                  selected from the pangenome for the 
                  similarity and also will prepare the 
                  blast files for the analysis")
p <- add_argument(p, "xgboost classifier", 
                  help = "if you have environmental
                  variables then you can train a classifier
                  with in the MetaPanRan using the XGBOOST
                  classification, both the baggin and the 
                  boosting based approaches have been 
                  implemented")
argv <- parse_args(p)
#####################################################################
############ Directory and dependencies check for processing ########
#####################################################################
if (!dir.exists("./roary_processing")) {
 dir.create("./roary_processing")
} 
if (!dir.exists("./peppan_processing")) {
 dir.create("./peppan_processing")
} 
if (!dir.exists("./pirate_processing")) {
 dir.create("./pirate_processing")
}
if (!dir.exists("./panroo_processing")) {
 dir.create("./panroo_processing")
}
if (!dir.exists("./prokka_annotations")) {
 dir.create("./prokka_annotations")
}
install_dep <- c("cdhit", "bamtools", "mafft", "fasttree", "diamond", 
                 "blast",  "iqtree", "muscle")
for (i in seq_along(install_dep)) {
 install_dep[i] <- system2("install_dep[i]", stdout = TRUE, stderr = TRUE)
} if (is.null(length(install_dep[i])) | file.access(install_dep[i]) == -1) {
 stop(paste0("dependencies for analysis are missing, please install dependencies
        using the following pattern 
        conda config --add channels defaults
        conda config --add channels bioconda
        conda config --add channels conda-forge
        conda install perl-bioperl==1.7.2 
        conda install mcl>=14.137 
        conda install mafft==7.310 
        conda install cd-hit>=4.6.4 
        conda install fasttree>=2.1.10 
        conda install diamond>=0.9.14 
        conda install blast>=2.2.31 
        conda install parallel>=20170422
        conda install mmseqs
        conda install rapidnj
        conda install iqtree
        conda install Fasttree
        brew install brewsci/bio/prokka
        brew install brewsci/bio/pirate"
 , sep = "\n"))
}
system2(command = "conda install", args = c("perl-bioperl", "mcl", "mafft", 
                       "cd-hit", "fasttree", "diamond", "blast","parallel", 
                              "mmseqs", "rapidnj", "iqtree", "Fasttree"))
system2(command =  "brew install brewsci/bio/pirate")
system2(command = "brew install brewsci/bio/prokka")
if (!file.exists("AMAS.py")) {
 system2(command = "pip install amas")
}
fasta_files <- list.files(path = "./prokka_annotations", recursive = TRUE, 
                          full.names = TRUE, pattern = "*.fasta")
for (i in seq_along(fasta_files)) {
 system2(command = "prokka", args = "fasta_files[i]", "--prefix fasta_files[i]")
 print("prokka annotations finished for the genomes and the path to the 
       annotations are in the directory prokka_annotations under the present
       working directory")
}
#################### roary run #####################################
if (argv$roary_run) {
 require(argv$gff_dir & argv$fasta_dir & 
          argv$metadata_path & argv$cluster_meta_path & argv$align)}
if (is.na(argv$gff_dir))
{stop("fasta directory containing the fasta files is missing")} 
if (is.na(argv$fasta_dir))
{stop("gff directory containing the gff files is missing")} 
if (is.na(argv$metadata_path))
{stop("meta data containing files is missing")} 
if (is.na(argv$cluster_meta_path))
{stop("meta data containing files is missing")}
gff_path <- argv$gff_dir
fasta_path <- argv$fasta_dir
meta_path <- argv$metadata_path
cluster_path <- argv$cluster_meta_path
align <- "argv$align"
system("cp -r gff_path/*.gff /roary_processing")
if (length(list.files(path = "roary_processing/",
              pattern = "*.gff" == length( list.files(path =
                                          "gff/path/", pattern = "*.gff")))))
{
 print("the gff files have been moved for the roary analysis
         and the file count has been verified")
}
roary_gff_files <- list.files( path = "./roary_processing",
                               pattern = "*.gff", recursive = TRUE, 
                               full.names = TRUE)
roary <- paste( "roary -e --mafft -i 90 -f ./roary_processing", 
                paste(roary_gff_files, collapse = " "))
system2(roary)
core_roary_alignment <- msaConvert(readDNAMultipleAlignment(
 filepath = "./roary_processing/core_gene_alignment.aln"), 
 type = "ape::DNAbin")
write.FASTA(core_roary_alignment, 
            file = "./roary_processing/core_roary_alignment.fasta")
core_roary_phylogeny <- system2(command = "iqtree", 
                       args = c("-s","core_roary_alignment.fasta" ,
                      "-p", "core_roary_alignment.phy" ,"-m" , "MFP+MERGE",
                                 "-B" , "1000", sep = " "))
core_roary_removed_gaps <- RemoveGaps(StaggerAlignment(AlignSeqs(
 readDNAStringSet(file = "./roary_processing/core_roary_alignment.fasta"))),
 removeGaps = "all", processors = 3)
core_roary_maximum_likelihood <- 
 TreeLine(myXStringSet = readDNAStringSet(file = 
                        "./roary_processing/core_roary_alignment.fasta"), 
        method = "ML", cutoff = 0.05, showPlot = TRUE, reconstruct = TRUE)
WriteDendrogram(core_roary_maximum_likelihood, 
                file = "./roary_processing/core_roary_maximum_likelihood.txt")
core_roary_genes_complete_likelihood <- TreeLine(myXStringSet = 
          readDNAStringSet(file = "core_roary_alignment.fasta"),myDistMatrix = 
         DistanceMatrix(RemoveGaps(StaggerAlignment(AlignSeqs(readDNAStringSet(
                     file = "./roary_processing/core_roary_alignment.fasta"))), 
                           removeGaps = "all", processors = 3), type = "dist"), 
       method = "complete", cutoff = 0.05, showPlot = TRUE, reconstruct = TRUE)
WriteDendrogram(core_roary_genes_complete_likelihood, file = "./roary_processing/core_roary_genes_complete_likelihood.txt")
alignment_plot <- ggmsa(DNAMultipleAlignment(msaMuscle(readDNAStringSet(
 file = "./roary_processing/core_roary_alignment.fasta"))), 
 color = "Shapely_NT", font = "DroidSansMono", 
 char_width = 0.5,seq_name = TRUE + geom_seqlogo(color = "Shapely_NT")) 
+geom_msaBar() + geom_GC()
ggsave(alignment_plot, width = 4, height = 4)
sink(file = "distance_matrix_of_the_roary_core_sequences.txt")
dist.dna(read.FASTA(file = "./roary_processing/core_genes_roary.fasta", 
                    type = "DNA"))
sink()
sink(file = "estimation_of_the_transition_and_transerversions_rate_of_the
            core_sequences.txt")
dnds(del.colgapsonly(msaMuscle(readDNAStringSet(
 file = "./roary_processing/core_roary_alignment.fasta"))))
sink()
roary_presence_absence <- list.files( path = "./roary_processing",
                         pattern = "*.csv", recursive = TRUE, full.names = TRUE)
roary_dataframe <- read.table(roary_presence_absence, header = TRUE, sep = "\t")
roary_fasta <- list.files(path = "fasta_path", recursive = TRUE, 
                          full.names = TRUE, pattern = "*.fasta")
fasta_final <- sapply(roary_fasta, readDNAStringSet)
roary_meta_read <- list.files( path = "meta_path", recursive = TRUE, 
                               full.names = TRUE, pattern = "*.meta.*.csv")
roary_meta <- read.table(roary_meta_read, header = TRUE, sep = "\t")
roary_cluster_read <- list.files( path = "cluster_path", recursive = TRUE,
                                  full.names = TRUE, pattern = "*.cluster.csv")
roary_cluster  <- read.table(roary_cluster_read, header = TRUE, sep = "\t")
roary_pagoo <- pagoo( data = roary_dataframe, org_meta = "roary_meta",
                      cluster_meta = "roary_cluster", sequences = "fasta_final", 
                      core_level = 100, sep = "__")
panviz_roary <- panviz(roary_presence_absence, location = "./roary_processing")

sequences <- roary_pagoo[["align"]]
write.FASTA(sequences, file = "./roary_processing/roary_selected_genes.fasta",
            header = TRUE)
roary_selected_genes_alignment <- msaConvert(msaMuscle(readDNAStringSet(
  file = "./roary_processing/roary_selected_genes.fasta"), type = "ape::DNAbin"))
write.FASTA(roary_selected_genes_alignment, 
            file = "./roary_processing/roary_selected_genes_alignment.fasta")                              
roary_selected_genes_alignment_consensus <- 
 msaConsensusSequence(msaMuscle(readDNAStringSet(
  file = "./roary_processing/roary_selected_genes.fasta")))
sink("./roary_processing/roary_selected_genes_alignment_consensus.txt")
msaConsensusSequence(msaMuscle(readDNAStringSet(
 file = "./roary_processing/roary_selected_genes.fasta")))
sink()
roary_selected_genes_alignment_consensus_data_frame <- 
 as.data.frame(consensusMatrix(msaConsensusSequence(msaMuscle(
  readDNAStringSet(file = "./roary_processing/roary_selected_genes.fasta")))),
  thresh = c(50, 20), ignoreGaps = TRUE)
sink(file = "roary_selected_genes_alignment_consensus_data_frame.txt")
as.character(consensusMatrix(msaConsensusSequence(msaMuscle(
 readDNAStringSet(file = "./roary_processing/roary_selected_genes.fasta")))),
 thresh = c(50, 20), ignoreGaps = TRUE)
sink()
roary_selected_genes_alignment_outlier <- 
 odseq(msaMuscle(readDNAStringSet(file = 
                             "./roary_processing/roary_selected_genes.fasta")), 
       distance_metric = "affine", B = 1000, threshold = 0.025)
sink(file = "roary_selected_genes_alignment_outlier.txt")
odseq(msaMuscle(readDNAStringSet(file = 
                             "./roary_processing/roary_selected_genes.fasta")), 
      distance_metric = "affine", B = 1000, threshold = 0.025)
sink()
core_sequences <- roary_pagoo$core_sequences_4_phylogeny(max_per_org = 1, 
                                                         fill = TRUE )
write.FASTA(core_sequences, file = 
             "./roary_processing/core_sequences_roary.fasta", header = TRUE)
read_core_sequences <- read.FASTA(file = 
                 "./roary_processing/core_sequences_roary.fasta", type = "DNA")
roary_core_sequences_alignment <- 
 msaConvert(msaMuscle(readDNAStringSet(file = 
        "./roary_processing/core_sequences_roary.fasta"),type = "ape::DNAbin"))
write.FASTA(roary_core_sequences_alignment, 
            file = "./roary_processing/roary_core_sequences_alignment.fasta")
roary_core_sequences_alignment_consensus <- 
 msaConsensusSequence(msaMuscle(readDNAStringSet(
  file = "./roary_processing/core_sequences_roary.fasta")))
sink(file = "roary_core_sequences_alignment_consensus.txt")
msaConsensusSequence(msaMuscle(readDNAStringSet(
 file = "./roary_processing/core_sequences_roary.fasta")))
sink()
roary_core_sequences_alignment_consensus_data_frame <- 
 as.data.frame(consensusMatrix(msaConsensusSequence(msaMuscle(
  readDNAStringSet(file = "./roary_processing/core_sequences_roary.fasta")))),
  thresh = c(50, 20), ignoreGaps = TRUE)
sink(file = "roary_core_sequences_alignment_consensus_data_frame.txt")
as.character(consensusMatrix(msaConsensusSequence(msaMuscle(
 readDNAStringSet(file = "./roary_processing/core_sequences_roary.fasta")))),
 thresh = c(50, 20), ignoreGaps = TRUE)
sink()
roary_core_sequences_phylogeny <- 
 system2(command = "iqtree", 
          args = c("-s","core_gene_roary_aligned.fasta",
                    "-p", "core_gene_roary_genes.phy", 
                             "-m" , "MFP+MERGE", "-B" , "1000", sep = " "))
roary_selected_genes_alignment_outlier <- odseq(msaMuscle(readDNAStringSet(
  file = "./roary_processing/core_gene_roary_aligned.fasta")),
  distance_metric = "affine", B = 1000, threshold = 0.025)
sink(file = "roary_selected_genes_alignment_outlier.txt")
odseq(msaMuscle(readDNAStringSet(file = 
                                  "./roary_processing/core_gene_roary_aligned.fasta")), 
      distance_metric = "affine", B = 1000, threshold = 0.025)
sink()
roary_core_genes_staggered_alignment <- 
 StaggerAlignment(AlignSeqs(readDNAStringSet(file = 
                                 "./roary_processing/core_genes_roary.fasta")))
roary_core_genes_staggered_alignment_write <- 
 msaConcvert(DNAMultipleAlignment(roary_core_genes_staggered_alignment), 
             type = "ape::DNAbin")
write.FASTA(roary_core_genes_staggered_alignment_write, 
            file = "./roary_processing/roary_core_genes_staggered_alignment_write.fasta")
roary_core_genes_staggered_alignment_removed_gaps <- RemoveGaps(StaggerAlignment(AlignSeqs(
 readDNAStringSet(file = "./roary_processing/core_genes_roary.fasta"))),
 removeGaps = "all", processors = 3)
roary_core_genes_staggered_alignment_removed_gaps_write <- 
 msaConcvert(DNAMultipleAlignment(roary_core_genes_staggered_alignment_removed_gaps), 
             type = "ape::DNAbin")
write.FASTA(roary_core_genes_staggered_alignment_write, 
  file = "./roary_processing/roary_core_genes_staggered_alignment_removed_gaps_write.fasta")
roary_core_genes_staggered_alignment_tidy <- tidy_msa(
 StaggerAlignment(AlignSeqs(readDNAStringSet(file = 
                                 "./roary_processing/core_genes_roary.fasta"))))
roary_core_genes_staggered_alignment_tidy_write <- 
 msaConcvert(DNAMultipleAlignment(roary_core_genes_staggered_alignment_tidy), 
             type = "ape::DNAbin")
write.FASTA(roary_core_genes_staggered_alignment_tidy_write, 
            file = "./roary_processing/roary_core_genes_staggered_alignment_tidy_write.fasta")
roary_core_genes_maximum_likelihood <- TreeLine(myXStringSet = 
          readDNAStringSet(file = "./roary_processing/core_genes_roary.fasta"), 
            method = "ML", cutoff = 0.05, showPlot = TRUE, reconstruct = TRUE)
WriteDendrogram(roary_core_genes_maximum_likelihood, 
                file = "./roary_processing/roary_core_genes_maximum_likelihood.txt")
roary_core_genes_complete_likelihood <- TreeLine(myXStringSet = 
              readDNAStringSet(file = "core_genes_roary.fasta"),myDistMatrix = 
      DistanceMatrix(RemoveGaps(StaggerAlignment(AlignSeqs(readDNAStringSet(
      file = "./roary_processing/core_genes_roary.fasta"))), removeGaps = "all", 
      processors = 3), type = "dist"), method = "complete", cutoff = 0.05, 
      showPlot = TRUE, reconstruct = TRUE)
WriteDendrogram(roary_core_genes_complete_likelihood, 
           file = "./roary_processing/roary_core_genes_complete_likelihood.txt")
alignment_plot <- 
 ggmsa(DNAMultipleAlignment(msaMuscle(readDNAStringSet(
  file = "core_genes_roary.fasta"))), color = "Shapely_NT", 
  font = "DroidSansMono", char_width = 0.5,seq_name = TRUE + 
   geom_seqlogo(color = "Shapely_NT")) +geom_msaBar() + geom_GC()
ggsave(alignment_plot, width = 4, height = 4)
pagoo_items <-   c( "pan_matrix", "organism", "clusters", "genes",
                    "sequences","core_level", "core_genes", "core_clusters",
                    "core_sequences", "cloud_genes","cloud_clusters",
                    "cloud_sequences", "shell_genes", "shell_clusters",
                    "summary_stats")
for (i in seq_along(pagoo_items))
{
 write.csv(roary_pagoo$pagoo_items[i], file =
            "roary_pagoo$pagoo_items[i].csv", sep = ",")
}
roary_pagoo$gg_pca(color = "organisms", size = 4) + theme_bw(base_size = 10)
+scale_color_brewer("Set2") + geom_point() + facet_grid("~")
distances = c("bray", "jaccard")
for (i in seq_along(distances)) {
 distances[i] <- roary_pagoo$dist(method = distances[i], binary = TRUE,
                                  diag = TRUE)
}
bar_plt <- roary_pagoo$gg_barplot()
bray_heatmap <- roary_pagoo$ggdist(method = "bray")
jaccard_heatmap <- roary_pagoo$ggdist(method = "jaccard")
roary_pangenome_curves <- roary_pagoo$gg_curves(what = "pangenome", size = 1) +
 ggtitle("pangenome curve of pirate_pagoo") +
 geom_point(alpha = 0.1, size = 4) +
 theme_classic() + ylim(0,50000) + scale_color_brewer("Accent")
roary_core_genome_curves <- roary_pagoo$gg_curves(what = "coregenome", size = 1) +
 ggtitle("pangenome curve of pirate_pagoo") +
 geom_point(alpha = 0.1, size = 4) +
 theme_classic() + ylim(0,50000) + scale_color_brewer(palette = "Set1") +
 scale_color_gradient() + scale_fill_brewer(direction = -1)
roary_pagoo$write_pangenome( dir = "roary_processing", force = FALSE)
ggsave("barplt.pdf", width = 4, height = 4)
ggsave("bray_heatmap",width = 4, height = 4 )
ggsave("jaccard_heatmap", width = 4, height = 4)
ggsave("roary_curves", width = 4, height = 4)
sink(file = "distance_matrix_of_the_roary_core_sequences.txt")
dist.dna(read.FASTA(file = "./roary_processing/core_genes_roary.fasta", 
                    type = "DNA"))
sink()
sink(file = "estimation_of_the_transition_and_transerversions_rate_of_the
            core_sequences.txt")
dnds(del.colgapsonly(msaMuscle(readDNAStringSet(
 file = "./roary_processing/core_genes_roary.fasta"))))
sink()
roary_selected_genes_for_kaks <- msaMuscle(
 readDNAStringSet(file =
              "./roary_processing/roary_selected_genes.fasta", format = "fasta")
)
roary_selected_aligned_gap_removed_genes_kaks <-
 msaConvert(RemoveGaps(DNAMultipleAlignment(roary_selected_genes_for_kaks)),
            type = "ape::DNAbin")
write.FASTA(roary_selected_aligned_gap_removed_genes_kaks,
      file = "./roary/processing/roary_selected_aligned_gap_removed_genes.fasta")
sink(file = "roary_selected_aligned_gap_removed_genes_kaks.txt")
kaks(read.alignment(file = "roary_selected_aligned_gap_removed_genes.fasta",
                format = "fasta"))
kaks(read.alignment(file = "roary_selected_aligned_gap_removed_genes.fasta",
                format = "fasta"))$ks
kaks(read.alignment(file = "roary_selected_aligned_gap_removed_genes.fasta",
                format = "fasta"))$ka
kaks(read.alignment(file = "roary_selected_aligned_gap_removed_genes.fasta",
                format = "fasta"))$vks
kaks(read.alignment(file = "roary_selected_aligned_gap_removed_genes.fasta",
                format = "fasta"))$vka
sink()
roary_core_genes_for_kaks <- msaMuscle(
 readDNAStringSet(file = "./roary_processing/core_sequences_roary.fasta",
                  format = "fasta")
)
roary_core_aligned_gap_removed_genes_kaks <- msaConvert(RemoveGaps(
 DNAMultipleAlignment(roary_core_genes_for_kaks)),type = "ape::DNAbin")
write.FASTA(roary_core_aligned_gap_removed_genes_kaks,
    file = "./roary/processing/roary_core_aligned_gap_removed_genes_kaks.fasta")
sink(file = "roary_selected_aligned_gap_removed_genes_kaks.txt")
kaks(read.alignment(file = "roary_core_aligned_gap_removed_genes_kaks.fasta",
                format = "fasta"))
kaks(read.alignment(file = "roary_core_aligned_gap_removed_genes_kaks.fasta",
                format = "fasta"))$ks
kaks(read.alignment(file = "roary_core_aligned_gap_removed_genes_kaks.fasta",
                format = "fasta"))$ka
kaks(read.alignment(file = "roary_core_aligned_gap_removed_genes_kaks.fasta",
                format = "fasta"))$vks
kaks(read.alignment(file = "roary_core_aligned_gap_removed_genes_kaks.fasta",
                format = "fasta"))$vka
sink()
roary_pan_phylogeny <- read.tree(TreeLine(myXStringSet =  
      readDNAStringSet(file = "./roary_processing/core_roary_alignment.fasta"), 
      method = "ML", cutoff = 0.05, showPlot = TRUE, reconstruct = TRUE))                         
roary_Rtab <- list.files(path = "./roary_processing", pattern = "*.Rtab")
roary_pan <- read.csv(roary_Rtab, row.names = 1, header = TRUE, sep = "\t")
pan_roary <- panstripe(roary_pan, roary_phylogeny)
sink(file = "write_panstripe_roary_summary.txt")
pan_roary$summary
sink()
plot_pangenome_params(pan_roary)
plot_pangenome_cumulative(pan_roary)
pan_roary_plot <- panstripe(roary_pan, roary_phylogeny, family = "gaussian")
plot_pangenome_params(panstr_gau) + plot_pangenome_cumulative(panstr_gau) + 
    plot_layout(nrows = 1)
roary_most_variable_genes <- colnames(roary_pan)[apply(roary_pan, 2, sd) > 6]
plot_tree_roary_variable_genes(tree = roary_phylogeny, pa = roary_pan, 
                    genes = variable_genes, label_genes = FALSE, col = "blacks")
plot_gain_loss(pan_roary)
plot_tsne(roary_pan)

roary_blast_read <- readDNAStringSet(file = 
                                "./roary/processing/core_genes_roary.fasta")
roary_sequences <- vector(mode = "character", length = length(roary_blast_read))

for (i in seq_along(roary_blast_read)){
  roary_sequences[i] <- roary_blast_read[i]
}

roary_blast_table <- list()
for (i in seq_along(rorary_sequences)){
  roary_blast_table[i] <- gget$blast(roary_sequence[i])
}
write.csv(roary_blast_table, file = "roary_blast_table.csv")
#################### PEPPAN run #####################################
if (argv$pEppan_run) {
 require(argv$gff_dir & argv$fasta_dir & 
          argv$metadata_path & argv$cluster_meta_path & argv$align)}
if (is.na(argv$gff_dir))
{stop("fasta directory containing the fasta files is missing")} 
if (is.na(argv$fasta_dir))
{stop("gff directory containing the gff files is missing")} 
if (is.na(argv$metadata_path))
{stop("meta data containing files is missing")} 
if (is.na(argv$cluster_meta_path))
{stop("meta data containing files is missing")}
gff_path <- argv$gff_dir
fasta_path <- argv$fasta_dir
meta_path <- argv$metadata_path
cluster_path <- argv$cluster_meta_path
align <- "argv$align"
system("cp -r gff_path/*.gff /peppan_processing")
if (length(list.files(path = "peppan_processing/",
                      pattern = "*.gff" == length( list.files(path =
                                   "gff/path/", pattern = "*.gff")))))
{
 print("the gff files have been moved for the peppan analysis
         and the file count has been verified")
}
peppan_gff_files <- list.files( path = "./peppan_processing",
                                pattern = "*.gff", recursive = TRUE, 
                                full.names = TRUE)
peppan <- paste( "PEPPAN", paste( "-p", peppan_gff_files[1], 
                                  roary_gff_files[-1], sep = " "))
system2(peppan)
system2("cd ./peppan_processing")
peppan_roary <- paste("python", "PEPPAN_parser.py - g","./peppan_processing/",
                      "*.PEPPA.gff", " - p", "peppan", " - m", " - t", " - a",
                      "100", " - ", "c", sep = " " )
system2(peppan_roary)
core_peppan_alignment <- msaConvert(readDNAMultipleAlignment(
 filepath = "./peppan_processing/core_gene_alignment.aln"), 
 type = "ape::DNAbin")
write.FASTA(core_peppan_alignment, 
            file = "./peppan_processing/core_peppan_alignment.fasta")
core_peppan_phylogeny <- system2(command = "iqtree", 
                                 args = c("-s","core_peppan_alignment.fasta" ,
                                          "-p", "core_peppan_alignment.phy" ,"-m" , "MFP+MERGE",
                                          "-B" , "1000", sep = " "))
core_peppan_removed_gaps <- RemoveGaps(StaggerAlignment(AlignSeqs(
 readDNAStringSet(file = "./peppan_processing/core_peppan_alignment.fasta"))),
 removeGaps = "all", processors = 3)
core_peppan_maximum_likelihood <- 
 TreeLine(myXStringSet = readDNAStringSet(file = 
                                           "./peppan_processing/core_peppan_alignment.fasta"), method = "ML", 
          cutoff = 0.05, showPlot = TRUE, reconstruct = TRUE)
WriteDendrogram(core_peppan_maximum_likelihood, 
                file = "./peppan_processing/core_peppan_maximum_likelihood.txt")
core_peppan_genes_complete_likelihood <- TreeLine(myXStringSet = 
          readDNAStringSet(file = "core_peppan_alignment.fasta"),myDistMatrix = 
         DistanceMatrix(RemoveGaps(StaggerAlignment(AlignSeqs(readDNAStringSet(
                  file = "./peppan_processing/core_peppan_alignment.fasta"))), 
                              removeGaps = "all", processors = 3), type = "dist"), 
        method = "complete", cutoff = 0.05, showPlot = TRUE, reconstruct = TRUE)
WriteDendrogram(core_peppan_genes_complete_likelihood, 
                file = "./peppan_processing/core_peppan_genes_complete_likelihood.txt")
alignment_plot <- ggmsa(DNAMultipleAlignment(msaMuscle(readDNAStringSet(
 file = "./peppan_processing/core_peppan_alignment.fasta"))), 
 color = "Shapely_NT", font = "DroidSansMono", 
 char_width = 0.5,seq_name = TRUE + geom_seqlogo(color = "Shapely_NT")) 
+geom_msaBar() + geom_GC()
ggsave(alignment_plot, width = 4, height = 4)
sink(file = "distance_matrix_of_the_peppan_core_sequences.txt")
dist.dna(read.FASTA(file = "./peppan_processing/core_genes_peppan.fasta", 
                    type = "DNA"))
sink()
sink(file = "estimation_of_the_transition_and_transerversions_rate_of_the
            core_sequences.txt")
dnds(del.colgapsonly(msaMuscle(readDNAStringSet(
 file = "./peppan_processing/core_peppan_alignment.fastaa"))))
sink()
peppan_presence_absence <- list.files( path = "./peppan_processing",
                                       pattern = "*.csv", recursive = TRUE, full.names = TRUE)
peppan_dataframe <- read.table(peppan_presence_absence, header = TRUE, sep = "\t")
peppan_fasta <- list.files(path = "fasta_path", recursive = TRUE, 
                           full.names = TRUE, pattern = "*.fasta")
fasta_final <- sapply(peppan_fasta, readDNAStringSet)
peppan_meta_read <- list.files( path = "meta_path", recursive = TRUE, 
                                full.names = TRUE, pattern = "*.meta.*.csv")
peppan_meta <- read.table(peppan_meta_read, header = TRUE, sep = "\t")
peppan_cluster_read <- list.files( path = "cluster_path", recursive = TRUE,
                                   full.names = TRUE, pattern = "*.cluster.csv")
peppan_cluster  <- read.table(peppan_cluster_read, header = TRUE, sep = "\t")
peppan_pagoo <- pagoo( data = peppan_dataframe, org_meta = "peppan_meta",
                       cluster_meta = "peppan_cluster", sequences = "fasta_final", 
                       core_level = 100, sep = "__")
panviz_peppan <- panviz(peppan_presence_absence, location = "./peppan_processing")

sequences <- peppan_pagoo[["align"]]
write.FASTA(sequences, file = "./peppan_processing/peppan_selected_genes.fasta",
            header = TRUE)
peppan_selected_genes_alignment <- msaConvert(msaMuscle(readDNAStringSet(
 file = "./peppan_processing/peppan_selected_genes.fasta"), 
                                             type = "ape::DNAbin"))
write.FASTA(peppan_selected_genes_alignment, 
            file = "./peppan_processing/peppan_selected_genes_alignment.fasta")                              
peppan_selected_genes_alignment_consensus <- 
 msaConsensusSequence(msaMuscle(readDNAStringSet(
  file = "./peppan_processing/peppan_selected_genes.fasta")))
sink("./peppan_processing/peppan_selected_genes_alignment_consensus.txt")
msaConsensusSequence(msaMuscle(readDNAStringSet(
 file = "./peppan_processing/peppan_selected_genes.fasta")))
sink()
peppan_selected_genes_alignment_consensus_data_frame <- 
 as.data.frame(consensusMatrix(msaConsensusSequence(msaMuscle(
  readDNAStringSet(file = "./peppan_processing/peppan_selected_genes.fasta")))),
  thresh = c(50, 20), ignoreGaps = TRUE)
sink(file = "peppan_selected_genes_alignment_consensus_data_frame.txt")
as.character(consensusMatrix(msaConsensusSequence(msaMuscle(
 readDNAStringSet(file = "./peppan_processing/peppan_selected_genes.fasta")))),
 thresh = c(50, 20), ignoreGaps = TRUE)
sink()
peppan_selected_genes_alignment_outlier <- 
 odseq(msaMuscle(readDNAStringSet(file = 
                                   "./peppan_processing/peppan_selected_genes.fasta")), 
       distance_metric = "affine", B = 1000, threshold = 0.025)
sink(file = "peppan_selected_genes_alignment_outlier.txt")
odseq(msaMuscle(readDNAStringSet(file = 
                                  "./peppan_processing/peppan_selected_genes.fasta")), 
      distance_metric = "affine", B = 1000, threshold = 0.025)
sink()
core_sequences <- peppan_pagoo$core_sequences_4_phylogeny(max_per_org = 1, 
                                                          fill = TRUE )
write.FASTA(core_sequences, file = 
             "./peppan_processing/core_sequences_peppan.fasta", header = TRUE)
read_core_sequences <- read.FASTA(file = 
                 "./peppan_processing/core_sequences_peppan.fasta", type = "DNA")
peppan_core_sequences_alignment <- 
 msaConvert(msaMuscle(readDNAStringSet(file = 
       "./peppan_processing/core_sequences_peppan.fasta"),type = "ape::DNAbin"))
write.FASTA(peppan_core_sequences_alignment, 
            file = "./peppan_processing/peppan_core_sequences_alignment.fasta")
peppan_core_sequences_alignment_consensus <- 
 msaConsensusSequence(msaMuscle(readDNAStringSet(
  file = "./peppan_processing/core_sequences_peppan.fasta")))
sink(file = "peppan_core_sequences_alignment_consensus.txt")
msaConsensusSequence(msaMuscle(readDNAStringSet(
 file = "./peppan_processing/core_sequences_peppan.fasta")))
sink()
peppan_core_sequences_alignment_consensus_data_frame <- 
 as.data.frame(consensusMatrix(msaConsensusSequence(msaMuscle(
  readDNAStringSet(file = "./peppan_processing/core_sequences_peppan.fasta")))),
  thresh = c(50, 20), ignoreGaps = TRUE)
sink(file = "peppan_core_sequences_alignment_consensus_data_frame.txt")
as.character(consensusMatrix(msaConsensusSequence(msaMuscle(
 readDNAStringSet(file = "./peppan_processing/core_sequences_peppan.fasta")))),
 thresh = c(50, 20), ignoreGaps = TRUE)
sink()
peppan_core_sequences_phylogeny <- 
 system2(command = "iqtree", args = c("-s","core_gene_peppan_aligned.fasta",
                                      "-p", "core_gene_peppan_genes.phy", 
                                      "-m" , "MFP+MERGE", "-B" , "1000", sep = " "))

peppan_selected_genes_alignment_outlier <- odseq(msaMuscle(readDNAStringSet(
 file = "./peppan_processing/core_gene_peppan_aligned.fasta")),
 distance_metric = "affine", B = 1000, threshold = 0.025) 
sink(file = "peppan_selected_genes_alignment_outlier.txt")
odseq(msaMuscle(readDNAStringSet(file = 
                                  "./peppan_processing/core_gene_peppan_aligned.fasta")), 
      distance_metric = "affine", B = 1000, threshold = 0.025)
sink()
peppan_core_genes_staggered_alignment <- 
 StaggerAlignment(AlignSeqs(readDNAStringSet(file = 
"./peppan_processing/core_genes_peppan.fasta")))
peppan_core_genes_staggered_alignment_write <- 
 msaConcvert(DNAMultipleAlignment(peppan_core_genes_staggered_alignment), 
             type = "ape::DNAbin")
write.FASTA(peppan_core_genes_staggered_alignment_write, 
            file = "./peppan_processing/peppan_core_genes_staggered_alignment_write.fasta")
peppan_core_genes_staggered_alignment_removed_gaps <- RemoveGaps(StaggerAlignment(AlignSeqs(
 readDNAStringSet(file = "./peppan_processing/core_genes_peppan.fasta"))),
 removeGaps = "all", processors = 3)
peppan_core_genes_staggered_alignment_removed_gaps_write <- 
 msaConcvert(DNAMultipleAlignment(peppan_core_genes_staggered_alignment_removed_gaps), 
type = "ape::DNAbin")
write.FASTA(peppan_core_genes_staggered_alignment_write, 
file = "./peppan_processing/peppan_core_genes_staggered_alignment_removed_gaps_write.fasta")
peppan_core_genes_staggered_alignment_tidy <- tidy_msa(
 StaggerAlignment(AlignSeqs(readDNAStringSet(file = 
"./peppan_processing/core_genes_peppan.fasta"))))
peppan_core_genes_staggered_alignment_tidy_write <- 
 msaConcvert(DNAMultipleAlignment(peppan_core_genes_staggered_alignment_tidy), 
             type = "ape::DNAbin")
write.FASTA(peppan_core_genes_staggered_alignment_tidy_write, 
     file = "./peppan_processing/peppan_core_genes_staggered_alignment_tidy_write.fasta")
peppan_core_genes_maximum_likelihood <- TreeLine(myXStringSet = 
readDNAStringSet(file = "./peppan_processing/core_genes_peppan.fasta"), 
method = "ML", cutoff = 0.05, showPlot = TRUE, reconstruct = TRUE)
WriteDendrogram(peppan_core_genes_maximum_likelihood, 
                file = "./peppan_processing/peppan_core_genes_maximum_likelihood.txt")
peppan_core_genes_complete_likelihood <- TreeLine(myXStringSet = 
readDNAStringSet(file = "core_genes_peppan.fasta"),myDistMatrix = 
DistanceMatrix(RemoveGaps(StaggerAlignment(AlignSeqs(readDNAStringSet(
file = "./peppan_processing/core_genes_peppan.fasta"))), removeGaps = "all", 
processors = 3), type = "dist"), method = "complete", cutoff = 0.05, 
showPlot = TRUE, reconstruct = TRUE)
WriteDendrogram(peppan_core_genes_complete_likelihood, 
                file = "./peppan_processing/peppan_core_genes_complete_likelihood.txt")
alignment_plot <- 
 ggmsa(DNAMultipleAlignment(msaMuscle(readDNAStringSet(
  file = "core_genes_peppan.fasta"))), color = "Shapely_NT", 
  font = "DroidSansMono", char_width = 0.5,seq_name = TRUE + 
   geom_seqlogo(color = "Shapely_NT")) +geom_msaBar() + geom_GC()
ggsave(alignment_plot, width = 4, height = 4)
pagoo_items <-   c( "pan_matrix", "organism", "clusters", "genes",
                    "sequences","core_level", "core_genes", "core_clusters",
                    "core_sequences", "cloud_genes","cloud_clusters",
                    "cloud_sequences", "shell_genes", "shell_clusters",
                    "summary_stats")
for (i in seq_along(pagoo_items))
{
 write.csv(peppan_pagoo$pagoo_items[i], file =
            "peppan_pagoo$pagoo_items[i].csv", sep = ",")
}
peppan_pagoo$gg_pca(color = "organisms", size = 4) + theme_bw(base_size = 10)
+scale_color_brewer("Set2") + geom_point() + facet_grid("~")
distances = c("bray", "jaccard")
for (i in seq_along(distances)) {
 distances[i] <- peppan_pagoo$dist(method = distances[i], binary = TRUE,
                                   diag = TRUE)
}
bar_plt <- peppan_pagoo$gg_barplot()
bray_heatmap <- peppan_pagoo$ggdist(method = "bray")
jaccard_heatmap <- peppan_pagoo$ggdist(method = "jaccard")
peppan_pangenome_curves <- peppan_pagoo$gg_curves(what = "pangenome", size = 1) +
 ggtitle("pangenome curve of pirate_pagoo") +
 geom_point(alpha = 0.1, size = 4) +
 theme_classic() + ylim(0,50000) + scale_color_brewer("Accent")
peppan_core_genome_curves <- peppan_pagoo$gg_curves(what = "coregenome", size = 1) +
 ggtitle("pangenome curve of pirate_pagoo") +
 geom_point(alpha = 0.1, size = 4) +
 theme_classic() + ylim(0,50000) + scale_color_brewer(palette = "Set1") +
 scale_color_gradient() + scale_fill_brewer(direction = -1)
peppan_pagoo$write_pangenome( dir = "peppan_processing", force = FALSE)
ggsave("barplt.pdf", width = 4, height = 4)
ggsave("bray_heatmap",width = 4, height = 4 )
ggsave("jaccard_heatmap", width = 4, height = 4)
ggsave("peppan_curves", width = 4, height = 4)
sink(file = "distance_matrix_of_the_peppan_core_sequences.txt")
dist.dna(read.FASTA(file = "./peppan_processing/core_genes_peppan.fasta", 
                    type = "DNA"))
sink()
sink(file = "estimation_of_the_transition_and_transerversions_rate_of_the
            core_sequences.txt")
dnds(del.colgapsonly(msaMuscle(readDNAStringSet(
 file = "./peppan_processing/core_genes_peppan.fasta"))))
sink()
peppan_selected_genes_for_kaks <- msaMuscle(
 readDNAStringSet(file =
"./peppan_processing/peppan_selected_genes.fasta", format = "fasta")
)
peppan_selected_aligned_gap_removed_genes_kaks <-
 msaConvert(RemoveGaps(DNAMultipleAlignment(peppan_selected_genes_for_kaks)),
 type = "ape::DNAbin")
write.FASTA(peppan_selected_aligned_gap_removed_genes_kaks,
 file = "./peppan/processing/peppan_selected_aligned_gap_removed_genes.fasta")
sink(file = "peppan_selected_aligned_gap_removed_genes_kaks.txt")
kaks(read.alignment(file = "peppan_selected_aligned_gap_removed_genes.fasta",
format = "fasta"))
kaks(read.alignment(file = "peppan_selected_aligned_gap_removed_genes.fasta",
  format = "fasta"))$ks
kaks(read.alignment(file = "peppan_selected_aligned_gap_removed_genes.fasta",
 format = "fasta"))$ka
kaks(read.alignment(file = "peppan_selected_aligned_gap_removed_genes.fasta",
format = "fasta"))$vks
kaks(read.alignment(file = "peppan_selected_aligned_gap_removed_genes.fasta",
format = "fasta"))$vka
sink()
peppan_core_genes_for_kaks <- msaMuscle(
 readDNAStringSet(file = "./peppan_processing/core_sequences_peppan.fasta",
format = "fasta")
)
peppan_core_aligned_gap_removed_genes_kaks <- msaConvert(RemoveGaps(
 DNAMultipleAlignment(peppan_core_genes_for_kaks)),type = "ape::DNAbin")

write.FASTA(peppan_core_aligned_gap_removed_genes_kaks,
file = "./peppan/processing/peppan_core_aligned_gap_removed_genes_kaks.fasta")
sink(file = "peppan_selected_aligned_gap_removed_genes_kaks.txt")
kaks(read.alignment(file = "peppan_core_aligned_gap_removed_genes_kaks.fasta",
format = "fasta"))
kaks(read.alignment(file = "peppan_core_aligned_gap_removed_genes_kaks.fasta",
format = "fasta"))$ks
kaks(read.alignment(file = "peppan_core_aligned_gap_removed_genes_kaks.fasta",
format = "fasta"))$ka
kaks(read.alignment(file = "peppan_core_aligned_gap_removed_genes_kaks.fasta",
format = "fasta"))$vks
kaks(read.alignment(file = "peppan_core_aligned_gap_removed_genes_kaks.fasta",
                    format = "fasta"))$vka
sink()
#################### PIRATE run #####################################
if (argv$pIrate_run) {
 require(argv$gff_dir & argv$fasta_dir & 
          argv$metadata_path & argv$cluster_meta_path & argv$align)}
if (is.na(argv$gff_dir))
{stop("fasta directory containing the fasta files is missing")} 
if (is.na(argv$fasta_dir))
{stop("gff directory containing the gff files is missing")} 
if (is.na(argv$metadata_path))
{stop("meta data containing files is missing")} 
if (is.na(argv$cluster_meta_path))
{stop("meta data containing files is missing")}
gff_path <- argv$gff_dir
fasta_path <- argv$fasta_dir
meta_path <- argv$metadata_path
cluster_path <- argv$cluster_meta_path
align <- "argv$align"
system("cp -r gff_path/*.gff /pirate_processing")
if (length(list.files(path = "pirate_processing/",
pattern = "*.gff" == length( list.files(path =
"gff/path/", pattern = "*.gff")))))
{
 print("the gff files have been moved for the pirate analysis
         and the file count has been verified")
}
pirate_gff_files <- list.files( path = "./pirate_processing",
                                pattern = "*.gff", recursive = TRUE, 
                                full.names = TRUE)
system2(command = "PIRATE", args = c("-i", "./pirate_processing", "-a", 
                "-r", "-k", "--diamond", "-s", "90,91,92,93,94,95", 
                "--cd-step 1", "--cd-low 95"))
system2(command = "PIRATE_to_roary.pl", args = 
         "-i", "./pirate_processing/PIRATE.*.tsv", 
        "-o", "./pirate_to_roary.tsv")
pirate_dataframe <- read.table("./pirate_processing/pirate_to_roary.tsv", 
                               header = TRUE, sep = "\t")
core_pirate_alignment <- msaConvert(readDNAMultipleAlignment(
 filepath = "./pirate_processing/core_gene_alignment.aln"), 
 type = "ape::DNAbin")
write.FASTA(core_pirate_alignment, 
            file = "./pirate_processing/core_pirate_alignment.fasta")
core_pirate_phylogeny <- system2(command = "iqtree", 
           args = c("-s","core_pirate_alignment.fasta" ,
            "-p", "core_pirate_alignment.phy" ,"-m" , "MFP+MERGE",
            "-B" , "1000", sep = " "))
core_pirate_removed_gaps <- RemoveGaps(StaggerAlignment(AlignSeqs(
 readDNAStringSet(file = "./pirate_processing/core_pirate_alignment.fasta"))),
 removeGaps = "all", processors = 3)
core_pirate_maximum_likelihood <- 
 TreeLine(myXStringSet = readDNAStringSet(file = 
             "./pirate_processing/core_pirate_alignment.fasta"), method = "ML", 
          cutoff = 0.05, showPlot = TRUE, reconstruct = TRUE)
WriteDendrogram(core_pirate_maximum_likelihood, 
                file = "./pirate_processing/core_pirate_maximum_likelihood.txt")
core_pirate_genes_complete_likelihood <- TreeLine(myXStringSet = 
       readDNAStringSet(file = "core_pirate_alignment.fasta"),myDistMatrix = 
       DistanceMatrix(RemoveGaps(StaggerAlignment(AlignSeqs(readDNAStringSet(
       file = "./pirate_processing/core_pirate_alignment.fasta"))), 
       removeGaps = "all", processors = 3), type = "dist"), 
       method = "complete", cutoff = 0.05, showPlot = TRUE, reconstruct = TRUE)
WriteDendrogram(core_pirate_genes_complete_likelihood, 
                file = "./pirate_processing/core_pirate_genes_complete_likelihood.txt")
alignment_plot <- ggmsa(DNAMultipleAlignment(msaMuscle(readDNAStringSet(
 file = "./pirate_processing/core_pirate_alignment.fasta"))), 
 color = "Shapely_NT", font = "DroidSansMono", 
 char_width = 0.5,seq_name = TRUE + geom_seqlogo(color = "Shapely_NT")) 
+geom_msaBar() + geom_GC()
ggsave(alignment_plot, width = 4, height = 4)
sink(file = "distance_matrix_of_the_pirate_core_sequences.txt")
dist.dna(read.FASTA(file = "./pirate_processing/core_genes_pirate.fasta", 
                    type = "DNA"))
sink()
sink(file = "estimation_of_the_transition_and_transerversions_rate_of_the
            core_sequences.txt")
dnds(del.colgapsonly(msaMuscle(readDNAStringSet(
 file = "./pirate_processing/core_pirate_alignment.fastaa"))))
sink()
pirate_presence_absence <- list.files( path = "./pirate_processing",
                                       pattern = "*.csv", recursive = TRUE, full.names = TRUE)
pirate_dataframe <- read.table(pirate_presence_absence, header = TRUE, sep = "\t")
pirate_fasta <- list.files(path = "fasta_path", recursive = TRUE, 
                           full.names = TRUE, pattern = "*.fasta")
fasta_final <- sapply(pirate_fasta, readDNAStringSet)
pirate_meta_read <- list.files( path = "meta_path", recursive = TRUE, 
                                full.names = TRUE, pattern = "*.meta.*.csv")
pirate_meta <- read.table(pirate_meta_read, header = TRUE, sep = "\t")
pirate_cluster_read <- list.files( path = "cluster_path", recursive = TRUE,
                                   full.names = TRUE, pattern = "*.cluster.csv")
pirate_cluster  <- read.table(pirate_cluster_read, header = TRUE, sep = "\t")
pirate_pagoo <- pagoo( data = pirate_dataframe, org_meta = "pirate_meta",
                       cluster_meta = "pirate_cluster", sequences = "fasta_final", 
                       core_level = 100, sep = "__")
panviz_pirate <- panviz(pirate_presence_absence, location = "./pirate_processing")

sequences <- pirate_pagoo[["align"]]
write.FASTA(sequences, file = "./pirate_processing/pirate_selected_genes.fasta",
            header = TRUE)
pirate_selected_genes_alignment <- msaConvert(msaMuscle(readDNAStringSet(file = 
             "./pirate_processing/pirate_selected_genes.fasta"), 
                                 type = "ape::DNAbin"))
write.FASTA(pirate_selected_genes_alignment, 
            file = "./pirate_processing/pirate_selected_genes_alignment.fasta")                              
pirate_selected_genes_alignment_consensus <- 
 msaConsensusSequence(msaMuscle(readDNAStringSet(
  file = "./pirate_processing/pirate_selected_genes.fasta")))
sink("./pirate_processing/pirate_selected_genes_alignment_consensus.txt")
msaConsensusSequence(msaMuscle(readDNAStringSet(
 file = "./pirate_processing/pirate_selected_genes.fasta")))
sink()
pirate_selected_genes_alignment_consensus_data_frame <- 
 as.data.frame(consensusMatrix(msaConsensusSequence(msaMuscle(
  readDNAStringSet(file = "./pirate_processing/pirate_selected_genes.fasta")))),
  thresh = c(50, 20), ignoreGaps = TRUE)
sink(file = "pirate_selected_genes_alignment_consensus_data_frame.txt")
as.character(consensusMatrix(msaConsensusSequence(msaMuscle(
 readDNAStringSet(file = "./pirate_processing/pirate_selected_genes.fasta")))),
 thresh = c(50, 20), ignoreGaps = TRUE)
sink()
pirate_selected_genes_alignment_outlier <- 
 odseq(msaMuscle(readDNAStringSet(file = 
               "./pirate_processing/pirate_selected_genes.fasta")), 
       distance_metric = "affine", B = 1000, threshold = 0.025)
sink(file = "pirate_selected_genes_alignment_outlier.txt")
odseq(msaMuscle(readDNAStringSet(file = 
         "./pirate_processing/pirate_selected_genes.fasta")), 
      distance_metric = "affine", B = 1000, threshold = 0.025)
sink()
core_sequences <- pirate_pagoo$core_sequences_4_phylogeny(max_per_org = 1, 
             fill = TRUE )
write.FASTA(core_sequences, file = 
             "./pirate_processing/core_sequences_pirate.fasta", header = TRUE)
read_core_sequences <- read.FASTA(file = 
        "./pirate_processing/core_sequences_pirate.fasta", type = "DNA")
pirate_core_sequences_alignment <- 
 msaConvert(msaMuscle(readDNAStringSet(file = 
      "./pirate_processing/core_sequences_pirate.fasta"),type = "ape::DNAbin"))
write.FASTA(pirate_core_sequences_alignment, 
            file = "./pirate_processing/pirate_core_sequences_alignment.fasta")
pirate_core_sequences_alignment_consensus <- 
 msaConsensusSequence(msaMuscle(readDNAStringSet(
  file = "./pirate_processing/core_sequences_pirate.fasta")))
sink(file = "pirate_core_sequences_alignment_consensus.txt")
msaConsensusSequence(msaMuscle(readDNAStringSet(
 file = "./pirate_processing/core_sequences_pirate.fasta")))
sink()
pirate_core_sequences_alignment_consensus_data_frame <- 
 as.data.frame(consensusMatrix(msaConsensusSequence(msaMuscle(
  readDNAStringSet(file = "./pirate_processing/core_sequences_pirate.fasta")))),
  thresh = c(50, 20), ignoreGaps = TRUE)
sink(file = "pirate_core_sequences_alignment_consensus_data_frame.txt")
as.character(consensusMatrix(msaConsensusSequence(msaMuscle(
 readDNAStringSet(file = "./pirate_processing/core_sequences_pirate.fasta")))),
 thresh = c(50, 20), ignoreGaps = TRUE)
sink()
pirate_core_sequences_phylogeny <- 
 system2(command = "iqtree", args = c("-s","core_gene_pirate_aligned.fasta",
        "-p", "core_gene_pirate_genes.phy", 
         "-m" , "MFP+MERGE", "-B" , "1000", sep = " "))
pirate_selected_genes_alignment_outlier <- odseq(msaMuscle(readDNAStringSet(file = 
         "./pirate_processing/core_gene_pirate_aligned.fasta")), 
           distance_metric = "affine", B = 1000, threshold = 0.025)
sink(file = "pirate_selected_genes_alignment_outlier.txt")
odseq(msaMuscle(readDNAStringSet(file = 
       "./pirate_processing/core_gene_pirate_aligned.fasta")), 
      distance_metric = "affine", B = 1000, threshold = 0.025)
sink()
pirate_core_genes_staggered_alignment <- 
 StaggerAlignment(AlignSeqs(readDNAStringSet(file = 
         "./pirate_processing/core_genes_pirate.fasta")))
pirate_core_genes_staggered_alignment_write <- 
 msaConcvert(DNAMultipleAlignment(pirate_core_genes_staggered_alignment), 
             type = "ape::DNAbin")
write.FASTA(pirate_core_genes_staggered_alignment_write, 
            file = "./pirate_processing/pirate_core_genes_staggered_alignment_write.fasta")
pirate_core_genes_staggered_alignment_removed_gaps <- RemoveGaps(StaggerAlignment(AlignSeqs(
 readDNAStringSet(file = "./pirate_processing/core_genes_pirate.fasta"))),
 removeGaps = "all", processors = 3)
pirate_core_genes_staggered_alignment_removed_gaps_write <- 
 msaConcvert(DNAMultipleAlignment(pirate_core_genes_staggered_alignment_removed_gaps), 
             type = "ape::DNAbin")
write.FASTA(pirate_core_genes_staggered_alignment_write, 
            file = "./pirate_processing/pirate_core_genes_staggered_alignment_removed_gaps_write.fasta")
pirate_core_genes_staggered_alignment_tidy <- tidy_msa(
 StaggerAlignment(AlignSeqs(readDNAStringSet(file = 
                                              "./pirate_processing/core_genes_pirate.fasta"))))
pirate_core_genes_staggered_alignment_tidy_write <- 
 msaConcvert(DNAMultipleAlignment(pirate_core_genes_staggered_alignment_tidy), 
             type = "ape::DNAbin")
write.FASTA(pirate_core_genes_staggered_alignment_tidy_write, 
            file = "./pirate_processing/pirate_core_genes_staggered_alignment_tidy_write.fasta")
pirate_core_genes_maximum_likelihood <- TreeLine(myXStringSet = 
          readDNAStringSet(file = "./pirate_processing/core_genes_pirate.fasta"), 
          method = "ML", cutoff = 0.05, showPlot = TRUE, reconstruct = TRUE)
WriteDendrogram(pirate_core_genes_maximum_likelihood, 
                file = "./pirate_processing/pirate_core_genes_maximum_likelihood.txt")
pirate_core_genes_complete_likelihood <- TreeLine(myXStringSet = 
       readDNAStringSet(file = "core_genes_pirate.fasta"),myDistMatrix = 
     DistanceMatrix(RemoveGaps(StaggerAlignment(AlignSeqs(readDNAStringSet(
     file = "./pirate_processing/core_genes_pirate.fasta"))), removeGaps = "all", 
     processors = 3), type = "dist"), method = "complete", cutoff = 0.05, 
     showPlot = TRUE, reconstruct = TRUE)
WriteDendrogram(pirate_core_genes_complete_likelihood, 
                file = "./pirate_processing/pirate_core_genes_complete_likelihood.txt")
alignment_plot <- 
 ggmsa(DNAMultipleAlignment(msaMuscle(readDNAStringSet(
  file = "core_genes_pirate.fasta"))), color = "Shapely_NT", 
  font = "DroidSansMono", char_width = 0.5,seq_name = TRUE + 
   geom_seqlogo(color = "Shapely_NT")) + geom_msaBar() + geom_GC()
ggsave(alignment_plot, width = 4, height = 4)
pagoo_items <-   c( "pan_matrix", "organism", "clusters", "genes",
                    "sequences","core_level", "core_genes", "core_clusters",
                    "core_sequences", "cloud_genes","cloud_clusters",
                    "cloud_sequences", "shell_genes", "shell_clusters",
                    "summary_stats")
for (i in seq_along(pagoo_items))
{
 write.csv(pirate_pagoo$pagoo_items[i], file =
            "pirate_pagoo$pagoo_items[i].csv", sep = ",")
}
pirate_pagoo$gg_pca(color = "organisms", size = 4) + theme_bw(base_size = 10)
+scale_color_brewer("Set2") + geom_point() + facet_grid("~")
distances = c("bray", "jaccard")
for (i in seq_along(distances)) {
 distances[i] <- pirate_pagoo$dist(method = distances[i], binary = TRUE,
                                   diag = TRUE)
}
bar_plt <- pirate_pagoo$gg_barplot()
bray_heatmap <- pirate_pagoo$ggdist(method = "bray")
jaccard_heatmap <- pirate_pagoo$ggdist(method = "jaccard")
pirate_pangenome_curves <- pirate_pagoo$gg_curves(what = "pangenome", size = 1) +
 ggtitle("pangenome curve of pirate_pagoo") +
 geom_point(alpha = 0.1, size = 4) +
 theme_classic() + ylim(0,50000) + scale_color_brewer("Accent")
pirate_core_genome_curves <- pirate_pagoo$gg_curves(what = "coregenome", size = 1) +
 ggtitle("pangenome curve of pirate_pagoo") +
 geom_point(alpha = 0.1, size = 4) +
 theme_classic() + ylim(0,50000) + scale_color_brewer(palette = "Set1") +
 scale_color_gradient() + scale_fill_brewer(direction = -1)
pirate_pagoo$write_pangenome( dir = "pirate_processing", force = FALSE)
ggsave("barplt.pdf", width = 4, height = 4)
ggsave("bray_heatmap",width = 4, height = 4 )
ggsave("jaccard_heatmap", width = 4, height = 4)
ggsave("pirate_curves", width = 4, height = 4)
sink(file = "distance_matrix_of_the_pirate_core_sequences.txt")
dist.dna(read.FASTA(file = "./pirate_processing/core_genes_pirate.fasta", 
                    type = "DNA"))
sink()
sink(file = "estimation_of_the_transition_and_transerversions_rate_of_the
            core_sequences.txt")
dnds(del.colgapsonly(msaMuscle(readDNAStringSet(
 file = "./pirate_processing/core_genes_pirate.fasta"))))
sink()
pirate_selected_genes_for_kaks <- msaMuscle(
 readDNAStringSet(file =
                   "./pirate_processing/pirate_selected_genes.fasta", format = "fasta")
)
pirate_selected_aligned_gap_removed_genes_kaks <-
 msaConvert(RemoveGaps(DNAMultipleAlignment(pirate_selected_genes_for_kaks)),
            type = "ape::DNAbin")
write.FASTA(pirate_selected_aligned_gap_removed_genes_kaks,
            file = "./pirate/processing/pirate_selected_aligned_gap_removed_genes.fasta")
sink(file = "pirate_selected_aligned_gap_removed_genes_kaks.txt")
kaks(read.alignment(file = "pirate_selected_aligned_gap_removed_genes.fasta",
                format = "fasta"))
kaks(read.alignment(file = "pirate_selected_aligned_gap_removed_genes.fasta",
                format = "fasta"))$ks
kaks(read.alignment(file = "pirate_selected_aligned_gap_removed_genes.fasta",
                format = "fasta"))$ka
kaks(read.alignment(file = "pirate_selected_aligned_gap_removed_genes.fasta",
                format = "fasta"))$vks
kaks(read.alignment(file = "pirate_selected_aligned_gap_removed_genes.fasta",
                format = "fasta"))$vka
sink()
pirate_core_genes_for_kaks <- msaMuscle(
 readDNAStringSet(file = "./pirate_processing/core_sequences_pirate.fasta",
                  format = "fasta")
)
pirate_core_aligned_gap_removed_genes_kaks <- msaConvert(RemoveGaps(
 DNAMultipleAlignment(pirate_core_genes_for_kaks)),type = "ape::DNAbin")

write.FASTA(pirate_core_aligned_gap_removed_genes_kaks,
            file = "./pirate/processing/pirate_core_aligned_gap_removed_genes_kaks.fasta")
sink(file = "pirate_selected_aligned_gap_removed_genes_kaks.txt")
kaks(read.alignment(file = "pirate_core_aligned_gap_removed_genes_kaks.fasta",
                format = "fasta"))
kaks(read.alignment(file = "pirate_core_aligned_gap_removed_genes_kaks.fasta",
                format = "fasta"))$ks
kaks(read.alignment(file = "pirate_core_aligned_gap_removed_genes_kaks.fasta",
                format = "fasta"))$ka
kaks(read.alignment(file = "pirate_core_aligned_gap_removed_genes_kaks.fasta",
                format = "fasta"))$vks
kaks(read.alignment(file = "pirate_core_aligned_gap_removed_genes_kaks.fasta",
                format = "fasta"))$vka
sink()
########################## Panroo run ###############################
if (argv$pAnroo_run) {
 require(argv$gff_dir & argv$fasta_dir & 
          argv$metadata_path & argv$cluster_meta_path & argv$align)}
if (is.na(argv$gff_dir))
{stop("fasta directory containing the fasta files is missing")} 
if (is.na(argv$fasta_dir))
{stop("gff directory containing the gff files is missing")} 
if (is.na(argv$metadata_path))
{stop("meta data containing files is missing")} 
if (is.na(argv$cluster_meta_path))
{stop("meta data containing files is missing")}
gff_path <- argv$gff_dir
fasta_path <- argv$fasta_dir
meta_path <- argv$metadata_path
cluster_path <- argv$cluster_meta_path
align <- "argv$align"
system("cp -r gff_path/*.gff /panroo_processing")
if (length(list.files(path = "panroo_processing/",
    pattern = "*.gff" == length( list.files(path =
              "gff/path/", pattern = "*.gff")))))
{
 print("the gff files have been moved for the panroo analysis
         and the file count has been verified")
}
panroo_gff_files <- list.files( path = "./panroo_processing",
                                pattern = "*.gff", recursive = TRUE, 
                                full.names = TRUE)
system2(command = "panaroo", args = ("*.gff", "-o", "results", 
                     "--clean-mode-strict", "--remove-invalid-genes"))
core_panroo_alignment <- msaConvert(readDNAMultipleAlignment(
 filepath = "./panroo_processing/core_gene_alignment.aln"), type = "ape::DNAbin")
write.FASTA(core_panroo_alignment, 
            file = "./panroo_processing/core_panroo_alignment.fasta")
core_panroo_phylogeny <- system2(command = "iqtree", 
       args = c("-s","core_panroo_alignment.fasta" ,
                     "-p", "core_panroo_alignment.phy" ,"-m" , "MFP+MERGE",
               "-B" , "1000", sep = " "))
core_panroo_removed_gaps <- RemoveGaps(StaggerAlignment(AlignSeqs(
 readDNAStringSet(file = "./panroo_processing/core_panroo_alignment.fasta"))),
 removeGaps = "all", processors = 3)
core_panroo_maximum_likelihood <- 
 TreeLine(myXStringSet = readDNAStringSet(file = 
             "./panroo_processing/core_panroo_alignment.fasta"), method = "ML", 
          cutoff = 0.05, showPlot = TRUE, reconstruct = TRUE)
WriteDendrogram(core_panroo_maximum_likelihood, 
                file = "./panroo_processing/core_panroo_maximum_likelihood.txt")
core_panroo_genes_complete_likelihood <- TreeLine(myXStringSet = 
           readDNAStringSet(file = "core_panroo_alignment.fasta"),myDistMatrix = 
        DistanceMatrix(RemoveGaps(StaggerAlignment(AlignSeqs(readDNAStringSet(
       file = "./panroo_processing/core_panroo_alignment.fasta"))), 
       removeGaps = "all", processors = 3), type = "dist"), 
   method = "complete", cutoff = 0.05, showPlot = TRUE, reconstruct = TRUE)
WriteDendrogram(core_panroo_genes_complete_likelihood, 
           file = "./panroo_processing/core_panroo_genes_complete_likelihood.txt")
alignment_plot <- ggmsa(DNAMultipleAlignment(msaMuscle(readDNAStringSet(
 file = "./panroo_processing/core_panroo_alignment.fasta"))), 
 color = "Shapely_NT", font = "DroidSansMono", 
 char_width = 0.5,seq_name = TRUE + geom_seqlogo(color = "Shapely_NT")) 
+geom_msaBar() + geom_GC()
ggsave(alignment_plot, width = 4, height = 4)
sink(file = "distance_matrix_of_the_panroo_core_sequences.txt")
dist.dna(read.FASTA(file = "./panroo_processing/core_genes_panroo.fasta", 
                    type = "DNA"))
sink()
sink(file = "estimation_of_the_transition_and_transerversions_rate_of_the
            core_sequences.txt")
dnds(del.colgapsonly(msaMuscle(readDNAStringSet(
 file = "./panroo_processing/core_panroo_alignment.fastaa"))))
sink()
panroo_presence_absence <- list.files( path = "./panroo_processing",
                                       pattern = "*.csv", recursive = TRUE, full.names = TRUE)
panroo_dataframe <- read.table(panroo_presence_absence, header = TRUE, sep = "\t")
panroo_fasta <- list.files(path = "fasta_path", recursive = TRUE, 
                           full.names = TRUE, pattern = "*.fasta")
fasta_final <- sapply(panroo_fasta, readDNAStringSet)
panroo_meta_read <- list.files( path = "meta_path", recursive = TRUE, 
                                full.names = TRUE, pattern = "*.meta.*.csv")
panroo_meta <- read.table(panroo_meta_read, header = TRUE, sep = "\t")
panroo_cluster_read <- list.files( path = "cluster_path", recursive = TRUE,
                                   full.names = TRUE, pattern = "*.cluster.csv")
panroo_cluster  <- read.table(panroo_cluster_read, header = TRUE, sep = "\t")
panroo_pagoo <- pagoo( data = panroo_dataframe, org_meta = "panroo_meta",
                       cluster_meta = "panroo_cluster", sequences = "fasta_final", 
                       core_level = 100, sep = "__")
panviz_panroo <- panviz(panroo_presence_absence, location = "./panroo_processing")
sequences <- panroo_pagoo[["align"]]
write.FASTA(sequences, file = "./panroo_processing/panroo_selected_genes.fasta",
            header = TRUE)
panroo_selected_genes_alignment <- msaConvert(msaMuscle(readDNAStringSet(file = 
            "./panroo_processing/panroo_selected_genes.fasta"), 
          type = "ape::DNAbin"))
write.FASTA(panroo_selected_genes_alignment, 
            file = "./panroo_processing/panroo_selected_genes_alignment.fasta")                              
panroo_selected_genes_alignment_consensus <- 
 msaConsensusSequence(msaMuscle(readDNAStringSet(
  file = "./panroo_processing/panroo_selected_genes.fasta")))
sink("./panroo_processing/panroo_selected_genes_alignment_consensus.txt")
msaConsensusSequence(msaMuscle(readDNAStringSet(
 file = "./panroo_processing/panroo_selected_genes.fasta")))
sink()
panroo_selected_genes_alignment_consensus_data_frame <- 
 as.data.frame(consensusMatrix(msaConsensusSequence(msaMuscle(
  readDNAStringSet(file = "./panroo_processing/panroo_selected_genes.fasta")))),
  thresh = c(50, 20), ignoreGaps = TRUE)
sink(file = "panroo_selected_genes_alignment_consensus_data_frame.txt")
as.character(consensusMatrix(msaConsensusSequence(msaMuscle(
 readDNAStringSet(file = "./panroo_processing/panroo_selected_genes.fasta")))),
 thresh = c(50, 20), ignoreGaps = TRUE)
sink()
panroo_selected_genes_alignment_outlier <- 
 odseq(msaMuscle(readDNAStringSet(file = 
                                   "./panroo_processing/panroo_selected_genes.fasta")), 
       distance_metric = "affine", B = 1000, threshold = 0.025)
sink(file = "panroo_selected_genes_alignment_outlier.txt")
odseq(msaMuscle(readDNAStringSet(file = 
                                  "./panroo_processing/panroo_selected_genes.fasta")), 
      distance_metric = "affine", B = 1000, threshold = 0.025)
sink()
core_sequences <- panroo_pagoo$core_sequences_4_phylogeny(max_per_org = 1, 
                                                          fill = TRUE )
write.FASTA(core_sequences, file = 
             "./panroo_processing/core_sequences_panroo.fasta", header = TRUE)
read_core_sequences <- read.FASTA(file = 
                                   "./panroo_processing/core_sequences_panroo.fasta", type = "DNA")
panroo_core_sequences_alignment <- 
 msaConvert(msaMuscle(readDNAStringSet(file = 
                                        "./panroo_processing/core_sequences_panroo.fasta"),type = "ape::DNAbin"))
write.FASTA(panroo_core_sequences_alignment, 
            file = "./panroo_processing/panroo_core_sequences_alignment.fasta")
panroo_core_sequences_alignment_consensus <- 
 msaConsensusSequence(msaMuscle(readDNAStringSet(
  file = "./panroo_processing/core_sequences_panroo.fasta")))
sink(file = "panroo_core_sequences_alignment_consensus.txt")
msaConsensusSequence(msaMuscle(readDNAStringSet(
 file = "./panroo_processing/core_sequences_panroo.fasta")))
sink()
panroo_core_sequences_alignment_consensus_data_frame <- 
 as.data.frame(consensusMatrix(msaConsensusSequence(msaMuscle(
  readDNAStringSet(file = "./panroo_processing/core_sequences_panroo.fasta")))),
  thresh = c(50, 20), ignoreGaps = TRUE)
sink(file = "panroo_core_sequences_alignment_consensus_data_frame.txt")
as.character(consensusMatrix(msaConsensusSequence(msaMuscle(
 readDNAStringSet(file = "./panroo_processing/core_sequences_panroo.fasta")))),
 thresh = c(50, 20), ignoreGaps = TRUE)
sink()
panroo_core_sequences_phylogeny <- 
 system2(command = "iqtree", args = c("-s","core_gene_panroo_aligned.fasta",
                                      "-p", "core_gene_panroo_genes.phy", 
                                      "-m" , "MFP+MERGE", "-B" , "1000", sep = " "))

panroo_selected_genes_alignment_outlier <- odseq(msaMuscle(readDNAStringSet(file = 
           "./panroo_processing/core_gene_panroo_aligned.fasta")), 
          distance_metric = "affine", B = 1000, threshold = 0.025)
sink(file = "panroo_selected_genes_alignment_outlier.txt")
odseq(msaMuscle(readDNAStringSet(file = 
                                  "./panroo_processing/core_gene_panroo_aligned.fasta")), 
      distance_metric = "affine", B = 1000, threshold = 0.025)
sink()
panroo_core_genes_staggered_alignment <- 
 StaggerAlignment(AlignSeqs(readDNAStringSet(file = 
                                              "./panroo_processing/core_genes_panroo.fasta")))
panroo_core_genes_staggered_alignment_write <- 
 msaConcvert(DNAMultipleAlignment(panroo_core_genes_staggered_alignment), 
             type = "ape::DNAbin")
write.FASTA(panroo_core_genes_staggered_alignment_write, 
            file = "./panroo_processing/panroo_core_genes_staggered_alignment_write.fasta")
panroo_core_genes_staggered_alignment_removed_gaps <- RemoveGaps(StaggerAlignment(AlignSeqs(
 readDNAStringSet(file = "./panroo_processing/core_genes_panroo.fasta"))),
 removeGaps = "all", processors = 3)
panroo_core_genes_staggered_alignment_removed_gaps_write <- 
 msaConcvert(DNAMultipleAlignment(panroo_core_genes_staggered_alignment_removed_gaps), 
             type = "ape::DNAbin")
write.FASTA(panroo_core_genes_staggered_alignment_write, 
            file = "./panroo_processing/panroo_core_genes_staggered_alignment_removed_gaps_write.fasta")
panroo_core_genes_staggered_alignment_tidy <- tidy_msa(
 StaggerAlignment(AlignSeqs(readDNAStringSet(file = 
                                              "./panroo_processing/core_genes_panroo.fasta"))))
panroo_core_genes_staggered_alignment_tidy_write <- 
 msaConcvert(DNAMultipleAlignment(panroo_core_genes_staggered_alignment_tidy), 
             type = "ape::DNAbin")
write.FASTA(panroo_core_genes_staggered_alignment_tidy_write, 
            file = "./panroo_processing/panroo_core_genes_staggered_alignment_tidy_write.fasta")
panroo_core_genes_maximum_likelihood <- TreeLine(myXStringSet = 
      readDNAStringSet(file = "./panroo_processing/core_genes_panroo.fasta"), 
       method = "ML", cutoff = 0.05, showPlot = TRUE, reconstruct = TRUE)
WriteDendrogram(panroo_core_genes_maximum_likelihood, 
                file = "./panroo_processing/panroo_core_genes_maximum_likelihood.txt")
panroo_core_genes_complete_likelihood <- TreeLine(myXStringSet = 
      readDNAStringSet(file = "core_genes_panroo.fasta"),myDistMatrix = 
      DistanceMatrix(RemoveGaps(StaggerAlignment(AlignSeqs(readDNAStringSet(
      file = "./panroo_processing/core_genes_panroo.fasta"))), removeGaps = "all", 
     processors = 3), type = "dist"), method = "complete", cutoff = 0.05, 
     showPlot = TRUE, reconstruct = TRUE)
WriteDendrogram(panroo_core_genes_complete_likelihood, 
                file = "./panroo_processing/panroo_core_genes_complete_likelihood.txt")
alignment_plot <- 
 ggmsa(DNAMultipleAlignment(msaMuscle(readDNAStringSet(
  file = "core_genes_panroo.fasta"))), color = "Shapely_NT", 
  font = "DroidSansMono", char_width = 0.5,seq_name = TRUE + 
   geom_seqlogo(color = "Shapely_NT")) + geom_msaBar() + geom_GC()
ggsave(alignment_plot, width = 4, height = 4)
pagoo_items <-   c( "pan_matrix", "organism", "clusters", "genes",
                    "sequences","core_level", "core_genes", "core_clusters",
                    "core_sequences", "cloud_genes","cloud_clusters",
                    "cloud_sequences", "shell_genes", "shell_clusters",
                    "summary_stats")
for (i in seq_along(pagoo_items))
{
 write.csv(panroo_pagoo$pagoo_items[i], file =
            "panroo_pagoo$pagoo_items[i].csv", sep = ",")
}
panroo_pagoo$gg_pca(color = "organisms", size = 4) + theme_bw(base_size = 10)
+scale_color_brewer("Set2") + geom_point() + facet_grid("~")
distances = c("bray", "jaccard")
for (i in seq_along(distances)) {
 distances[i] <- panroo_pagoo$dist(method = distances[i], binary = TRUE,
                                   diag = TRUE)
}
bar_plt <- panroo_pagoo$gg_barplot()
bray_heatmap <- panroo_pagoo$ggdist(method = "bray")
jaccard_heatmap <- panroo_pagoo$ggdist(method = "jaccard")
panroo_pangenome_curves <- panroo_pagoo$gg_curves(what = "pangenome", size = 1) +
 ggtitle("pangenome curve of panroo_pagoo") +
 geom_point(alpha = 0.1, size = 4) +
 theme_classic() + ylim(0,50000) + scale_color_brewer("Accent")
panroo_core_genome_curves <- panroo_pagoo$gg_curves(what = "coregenome", size = 1) +
 ggtitle("pangenome curve of panroo_pagoo") +
 geom_point(alpha = 0.1, size = 4) +
 theme_classic() + ylim(0,50000) + scale_color_brewer(palette = "Set1") +
 scale_color_gradient() + scale_fill_brewer(direction = -1)
panroo_pagoo$write_pangenome( dir = "panroo_processing", force = FALSE)
ggsave("barplt.pdf", width = 4, height = 4)
ggsave("bray_heatmap",width = 4, height = 4 )
ggsave("jaccard_heatmap", width = 4, height = 4)
ggsave("panroo_curves", width = 4, height = 4)
sink(file = "distance_matrix_of_the_panroo_core_sequences.txt")
dist.dna(read.FASTA(file = "./panroo_processing/core_genes_panroo.fasta", 
                    type = "DNA"))
sink()
sink(file = "estimation_of_the_transition_and_transerversions_rate_of_the
            core_sequences.txt")
dnds(del.colgapsonly(msaMuscle(readDNAStringSet(
 file = "./panroo_processing/core_genes_panroo.fasta"))))
sink()
panroo_selected_genes_for_kaks <- msaMuscle(
 readDNAStringSet(file =
                   "./panroo_processing/panroo_selected_genes.fasta", format = "fasta")
)
panroo_selected_aligned_gap_removed_genes_kaks <-
 msaConvert(RemoveGaps(DNAMultipleAlignment(panroo_selected_genes_for_kaks)),
            type = "ape::DNAbin")
write.FASTA(panroo_selected_aligned_gap_removed_genes_kaks,
            file = "./panroo/processing/panroo_selected_aligned_gap_removed_genes.fasta")
sink(file = "panroo_selected_aligned_gap_removed_genes_kaks.txt")
kaks(read.alignment(file = "panroo_selected_aligned_gap_removed_genes.fasta",
                format = "fasta"))
kaks(read.alignment(file = "panroo_selected_aligned_gap_removed_genes.fasta",
                format = "fasta"))$ks
kaks(read.alignment(file = "panroo_selected_aligned_gap_removed_genes.fasta",
                format = "fasta"))$ka
kaks(read.alignment(file = "panroo_selected_aligned_gap_removed_genes.fasta",
                format = "fasta"))$vks
kaks(read.alignment(file = "panroo_selected_aligned_gap_removed_genes.fasta",
                format = "fasta"))$vka
sink()
panroo_core_genes_for_kaks <- msaMuscle(
 readDNAStringSet(file = "./panroo_processing/core_sequences_panroo.fasta",
                  format = "fasta")
)
panroo_core_aligned_gap_removed_genes_kaks <- msaConvert(RemoveGaps(
 DNAMultipleAlignment(panroo_core_genes_for_kaks)),type = "ape::DNAbin")
write.FASTA(panroo_core_aligned_gap_removed_genes_kaks,
            file = "./panroo/processing/panroo_core_aligned_gap_removed_genes_kaks.fasta")
sink(file = "panroo_selected_aligned_gap_removed_genes_kaks.txt")
kaks(
 read.alignment(file = "panroo_core_aligned_gap_removed_genes_kaks.fasta",
                format = "fasta"))
kaks(read.alignment(file = "panroo_core_aligned_gap_removed_genes_kaks.fasta",
                format = "fasta"))$ks
kaks(read.alignment(file = "panroo_core_aligned_gap_removed_genes_kaks.fasta",
                format = "fasta"))$ka
kaks(read.alignment(file = "panroo_core_aligned_gap_removed_genes_kaks.fasta",
                format = "fasta"))$vks
kaks(read.alignment(file = "panroo_core_aligned_gap_removed_genes_kaks.fasta",
                format = "fasta"))$vka
sink()
############ Venn visualization of the common metagenomes ###########
if (argv$venn_viz) {
 require("roary_venn", "peppan_venn", "pirate_venn", "panaroo_venn" )
 stop(" the venn visualiation requires roary_presence_absence,
        peppan_presence_absence, pirate_presence_absence, 
        panaroo_presence_absence")
} 
roary_venn <-  list.files( path = "./roary_processing", recursive = TRUE, 
                           full.names = TRUE, pattern = "*.roary.csv")
peppan_venn <- list.files( path = "./peppan_processing", recursive = TRUE,
                           full.names = TRUE, pattern = "*.peppan.csv")
pirate_venn <- list.files( path = "./pirate_processing", recursive = TRUE,
                           full.names = TRUE, pattern = "*.pirate.csv")
panroo_venn <- list.files( path = "./pirate_processing", recursive = TRUE,
                           full.names = TRUE, pattern = "*.panroo.csv")
roary_organism <- roary_venn$organisms
peppan_organism <- peppan_venn$organisms
pirate_organism <- pirate_venn$organisms
panroo_organism <- panroo_venn$organisms
matrix(roary_organism, peppan_organism, pirate_organism, panroo_organism)
venn_dataframe <- as.data.frame(matrix)
venn(venn_dataframe, col = "navyblue", ellipses = TRUE, 
     ilabels = TRUE, zcolor = "style")

#citation
# MetaPanRan: A standalone R based integrated analysis for metagenomes from 
# pangenomics to machine learning. 
# Gaurav Sablok, Thulani Makhalanyane
# Department of Biochemistry, Genetics and Microbiology,
# Faculty of Natural and Agricultural Sciences, 
# UNIVERSITY OF PRETORIA
# read manual
# Abstract: 
#   Owing to the wide diversity and growing enormous potential of metagenomics 
# widely associated with many of the biological spans, it is detrimental to have 
# potential look at the meta pangenomics. In the present paper, we present an 
# integrated work-flow for analyzing the microbial pangenomes using several state 
# of the art pangenomics methods such as ROARY, PANAROO, PIRATE and PEPPAN. 
# MetaPanRan allows you to run all these methods as a part of the application 
# and allows for the downstream analysis of the microbial pangenomes such as the 
# alignments, core alignments, alignment filtering, phylogeny, estimation of dn/ds, 
# detection of outliers, and visualization and estimation of pangenomes. 
# MetaPanRan serve as a common platform for the detection of the 
# metagenomics and downstream analysis associated with metagenomes. 
# 
# Microbes play an important role in shaping the diversity and classification of 
# the organisms as a result reflecting changes in the health as well as the diversity.
# Application of next generation sequencing methods has not only widely elaborated 
# the state of the art metagenomics by illustrating individual microbial strains as 
# well as the hierarchy level of the metagenomics at the strain level. Recently, a 
# surge of metagenomics application has seen in understanding the pangenomics of the 
# microbes, which has allowed for identifying common traits of interest as well as 
# parallel application of genomic breadth of approaches at a multi-level class and 
# at the phylum level. Several state of the art tools have been developed for the 
# identification of the metagenomes such as ROARY (Reference), PANAROO (Reference), 
# PIRATE (Reference), PEPPAN (Reference) and individual packages such as micropan 
# (Reference), Pagoo (Reference), Panstripe (Reference) have been developed for 
# the downstream analysis. However, the common limitation stays with the non-availability 
# of the streamline work-flow for analyzing the pangenomes all together.

# In this paper, we present MetagenomeIntegrated which addresses this problem by 
# providing a streamline work-flow by integrating all these metagenomes prediction 
# tools. In addition, MetaPanRan also provide options for doing the genome predictions 
# using the Prokka(Reference) and can be combined with the additional provided user 
# supplied annotation files (format = “*.gff”) for the downstream pangenome 
# predictions and annotations. MetaPanRan work-flow is depicted in Figure 1. It allows 
# the users to predict the meta-pangenomes using the ROARY (Reference), PANAROO (Reference),
# PIRATE (Reference), PEPPAN (Reference). Following the prediction of the metagenomes, 
# it performs several downstream analysis such as doing the pangenome alignments, 
# alignment visualization, masking the alignment for the gaps, and implementing a model 
# based phylogeny on the selected core genes and the pangenome alignments. It uses 
# muscle for the alignment and also for the staggered alignment and IQTREE(Reference) 
# under the “MFP+MERGE” model with 1000 bootstraps for the clade phylogeny. Maximum 
# likelihood phylogeny has been the main stream of development of the phylogeny with 
# or without constraints. MetaPanRan implements the maximum likelihood based phylogeny 
# using TreeLine with a cutoff of 0.5 (Reference). 

# Sequence visualization and looking at the sequence plots with the base probabilities 
# is one of the main stream of looking at the consensus alignments. MetaPanRan implements 
# ggmsa (Reference) for the selection and visualization of the base probabilities. 
# Panstripe (Reference) has been integrated in MetaPanRan for the visualization of 
# the core pan phylogeny and plotting the “gaussian” distribution of the pan metagenomes. 
# MetaPanRan requires minimal installation and  it performs all the installation and 
# directory wise hierarchical classification as a part of the pan genome predictions 
# and downstream analysis. All the core and pan genes files as a part of the user selected 
# genes are written as fasta formatted files and base probabilities. To furthermake it 
# more meaningful to the analysis, we integrated an in-build BLAST analysis using the 
# package ggET(Reference) which will blast all of the core genes and prepare a BLAST 
# table for the analysis to further check the sequence homology. 

# Microbial genes shows accelerated rate of evolution and this has been widely 
# demonstrated and has been used for simulation and genetic algorithms development. 
# To understand the accelerated rate of evolution, MetagenomeIntegrated first estimates 
# the alignment using Muscle(Reference) and then delete cols with gaps using the APE 
# package(Reference) and subsequently detected the outliers in the alignment using 
# OdSeq(Reference). Following the outlier detection and alignment filtering, it estimates 
# the transition/transversion rates using the seqinr package(Reference) and write a 
# dataframe object. To conclude, we present an integrated version of metagenome predictions, 
# analysis and visualization as a stand alone R tool, which allows the users to do meta 
# pangenome analysis with minimal interactions and at a fine scale of the metagenome. 
# This bridges a lack of the integrated workflow for the meta pangenomes. 