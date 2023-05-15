!/usr/bin/env Rscript --vanilla

# #####################################################################
# ################### Usage of Metagenome Integrated #################
# ####################################################################
# 
# (base) ➜  CodeR git:(main) ✗ ./metagenomics_Integrated.R -h
# usage: metagenomics_pangenome.R [--] [--help] [--opts OPTS]
# [--annoatate_genome ANNOATATE_GENOME] [--roary_run ROARY_RUN]
# [--PEPPAN_run PEPPAN_RUN] [--PIRATE_run PIRATE_RUN]
# [--PANAROO_file PANAROO_FILE] [--gff_dir GFF_DIR]
# [--metadata_path METADATA_PATH] [--cluster_meta_path
#  CLUSTER_META_PATH] [--fasta_dir FASTA_DIR] [--export EXPORT]
# [--align ALIGN] [--core_genes_phylogeny CORE_GENES_PHYLOGENY]
# [--venn_visualization VENN_VISUALIZATION]
# 
# An integrated pipeline for the prediction of the bacterial genomes and
# analysis and visualization of the pangenomes from ROARY, PIRATE,
# PEPPAN, PANROO. This allows you to run the three accurate pangenome
# predictions from ROARY, PEPPAN, PIRATE, PANROO and then gives the
# pangenomes. It also gives a comparative assessment of the pangenomes
# that is conserved across all through a venn visualization
# 
# flags:
#  -h, --help                show this help message and exit
# 
# optional arguments:
#  -x, --opts                RDS file containing argument values
# -a, --annoatate_genome    If you have just assembled your genome your
#                           can call this option to annotate the genome
# -r, --roary_run           Run roary as a part of the analysis
# -P, --PEPPAN_run          Will run the PEPPAN analysis and will
#                           convert to the roary for pangenome analysis
# --PIRATE_run              Will run the PIRATE analysis and will
#                           Convert from PIRATE output to roary for
#                           pangenome analysis
# --PANAROO_file            Process the panraoo file
# -g, --gff_dir             PATH to the gff directory to look for the
#                           presenece and absence gene
# -m, --metadata_path       add metadata for the pangenome
# -c, --cluster_meta_path   path for the meta data for the clusters
# -f, --fasta_dir           Provide the path to the folder containing
#                           the fasta file. You can provide the path to
#                           the directory
# -e, --export              Export the alignment for the particular
#                           genes or complete dataframe
# --align                   align the selected genes using the IQTREE
#                           and the using the MFP+MERGE model
#                           by finding the best partition scheme
#                           followed by tree inference and bootstrap
#                           and will also do a outlier detection in the
#                           aligned sequences and will plot a phylogeny
# --core_genes_phylogeny    this options allows to build the phylogeny
#                           using the best partition scheme followed by
#                           tree inference and bootstrap using all the
#                           core genes
# -v, --venn_viz            create a venn visualiztion of the common
#                           pangenome across all the roary, pirate,
#                           peppan, panroo, For this you have to
#                           provide the presence and absence file from
#                           each of the run after running the entire
#                           workflow or you can provide directly the
#                            file.
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


p <- arg_parser(
 " An integrated pipeline for the prediction of the bacterial genomes and 
    analysis and visualization of the pangenomes from ROARY, PIRATE, PEPPAN, 
    PANROO. This allows you to run the three accurate pangenome predictions 
    from ROARY, PEPPAN, PIRATE, PANROO and then gives the pangenomes. It also 
    gives a comparative assessment of the pangenomes that is conserved across 
    all through a venn visualization"
)
p <- add_argument(p, "--annoatate_genome", 
                  help = "If you have just assembled your genome
                  your can call this option to annotate the genome")
p <- add_argument(p, "--roary_run",
                  help = "Run roary as a part of the analysis")
p <- add_argument(p, "--PEPPAN_run",
                  help = "Will run the PEPPAN analysis and
                       will convert to the roary for pangenome 
                       analysis")
p <- add_argument(p, "--PIRATE_run",
                  help = "Will run the PIRATE analysis and 
                  will Convert from PIRATE output to roary
                  for pangenome analysis")
p <- add_argument(p, "--PANAROO_file",
                  help = "Process the panraoo file")
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
p <- add_argument(p, "--export",
                  help = "Export the alignment for the
                  particular genes or complete dataframe")
p <- add_argument(p, "--align",
                  help = "align the selected genes using
                  the IQTREE and the fastTREE using the 
                  MFP+MERGE model by finding the best 
                  partition scheme followed by tree 
                  inference and bootstrap and will also
                  do a outlier detection in the aligned
                  sequences and will plot a phylogeny")
p <- add_argument(p, "--core_genes_phylogeny", 
                  help = "this options allows to build 
                  the phylogeny using the best partition 
                  scheme followed by tree inference and 
                  bootstrap using all the core genes")
p <- add_argument(p, "--venn_viz", 
                  help = "create a venn visualiztion of 
                  the common pangenome across all the 
                  roary, pirate, peppan, panroo, For this
                  you have to provide the presence and
                  absence file from each of the run
                  after running the entire workflow
                  or you can provide directly the file.")
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
 stop( "dependencies for analysis are missing, please install dependencies
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
       )
}
system2(command = "conda install", args = c("perl-bioperl", "mcl", "mafft", 
                   "cd-hit", "fasttree", "diamond", "blast","parallel", 
                   "mmseqs", "rapidnj", "iqtree", "Fasttree"))
system2(command =  "brew install brewsci/bio/pirate")
system2(command = "brew install brewsci/bio/prokka")
#####################################################################
#################### Prokka Annotations #############################
#####################################################################
fasta_files <- list.files(path = "./prokka_annotations", recursive = TRUE, 
                full.names = TRUE, pattern = "*.fasta")
for (i in seq_along(fasta_files)) {
 system2(command = "prokka", args = "fasta_files[i]", "--prefix fasta_files[i]")
 print("prokka annotations finished for the genomes and the path to the 
       annotations are in the directory prokka_annotations under the present
       working directory")
}
#####################################################################
#################### Roary run ######################################
#####################################################################
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
panviz(roary_presence_absence, location = "./roary_processing")
sequences <- roary_pagoo[["align"]]
write.FASTA(sequences, file = "./roary_processing/selected_genes.fasta",
                                            header = TRUE)
alignment <- read.FASTA(file = "./roary_processing/selected_genes.fasta",
                        type = "DNA")
alignment <- msa(alignment)
write.FASTA(alignment, file = "./roary_processing/selected_genes_aligned.fasta",
                     header = TRUE)
selected_phylogeny <- system2(command = "iqtree", args =
                              c("-s", "selected_genes_aligned.fasta" ,
                          "-p", "selected_genes.phy" ,"-m" , "MFP+MERGE",
                          "-B" , "1000", sep = " "))
outlier <- odseq(alignment, distance_metric = "affine", B = 1000,
                 threshold = 0.025)
alignment_plot <- ggmsa( alignment, start, end, color = "Shapely_NT",
                         font = "DroidSansMono", char_width = 0.5,
                         seq_name = TRUE) + geom_seqlogo(color = "Shapely_NT")
core_sequences <- roary_pagoo$core_sequences_4_phylogeny
                                     (max_per_org = 1, fill = TRUE )
write.FASTA(core_sequences, file = "core_genes_roary.fasta", header = TRUE)
read_core_genes <- read.FASTA(file = "core_genes_roary.fasta", type = "DNA")
core_alignment <- msa(read_core_genes)
write.FASTA(core_alignment, file = "core_gene_roary_aligned.fasta",
                    header = TRUE)
core_phylogeny <- system2(command = "iqtree", args = c("-s",
                                 "core_gene_roary_aligned.fasta" ,
                          "-p", "selected_genes.phy" ,"-m" , "MFP+MERGE",
                          "-B" , "1000", sep = " "))
core_outlier <- odseq(core_alignment, distance_metric = "affine", B = 1000,
                 threshold = 0.025)
core_alignment_plot <- ggmsa( core_alignment, start, end, color = "Shapely_NT",
                         font = "DroidSansMono", char_width = 0.5,
                         seq_name = TRUE) + geom_seqlogo(color = "Shapely_NT")
+geom_msaBar() + geom_GC()
save.image(file = "roary_alignment_plot", version = 1)
save.image(file = "core_roary_alignment_plot", version = 1)
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
bar_plot <- roary_pagoo$gg_barplot()
bray_heatmap <- roary_pagoo$ggdist(method = "bray")
jaccard_heatmap <- roary_pagoo$ggdist(method = "jaccard")
roary_pangenome_curves <- roary_pagoo$gg_curves(what = "pangenome", size = 1) +
 ggtitle("pangenome curve of pirate_pagoo") +
 geom_point(alpha = 0.1, size = 4) +
 theme_classic() + ylim(0,50000) + scale_color_brewer("Accent")
roary_coregenome_curves <- roary_pagoo$gg_curves(what = "coregenome", size = 1) +
 ggtitle("pangenome curve of pirate_pagoo") +
 geom_point(alpha = 0.1, size = 4) +
 theme_classic() + ylim(0,50000) + scale_color_brewer(palette = "Set1") +
    scale_color_gradient() + scale_fill_brewer(direction = -1)
roary_pagoo$write_pangenome( dir = "roary_processing", force = FALSE)
save.image(file = "bar_plt", version = 1)
save.image(file = "bray_heatmap", version = 1)
save.image(file = "jaccard_heatmap", version = 1)
save.image(file = "roary_curves", version = 1)
# #####################################################################
# #################### PEPPAN run #####################################
# #####################################################################
if (argv$PEPPAN_run) {
 require(argv$fasta_dir & argv$gff_dir & argv$metadata_path & 
          argv$cluster_meta_path)}
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
align <- argv$align
system("cp -r gff_path/*.gff /peppan_analysis")
if (length(list.files(path = "peppan_analysis/",
                      pattern = "*.gff" == length(list.files(path =
                                  "gff/path/", pattern = "*.gff")))))
{
 print("the gff files have been moved for the roary analysis
         and the file count has been verified")
}
peppan_gff_files <- list.files( path = "./peppan_processing",
                      pattern = "*.gff", recursive = TRUE, full.names = TRUE)
peppan <- paste( "PEPPAN", paste( "-p", peppan_gff_files[1], 
                                  roary_gff_files[-1], sep = " "))
system2(peppan)
system2("cd ./peppan_processing")
peppan_roary <- paste("python", "PEPPAN_parser.py - g" "./peppan_processing/",
                      "*.PEPPA.gff", " - p", "peppan", " - m", " - t", " - a",
                      "100", " - ", "c", sep = " " )
system2(peppan_roary)
system2("cd ../")
peppan_presence_absence <- list.files( path = "./peppan_processing",
                                       pattern = "*.csv", recursive = TRUE, 
                                       full.names = TRUE)
peppan_dataframe <- read.table(peppan_presence_absence, header = TRUE, 
                               sep = "\t")
peppan_fasta <- list.files(path = "fasta_path", recursive = TRUE, 
                           full.names = TRUE, pattern = "*.fasta")
fasta_final <- sapply(peppan_fasta, readDNAStringSet)
peppan_meta_read <- list.files( path = "meta_path", recursive = TRUE, 
                                full.names = TRUE, pattern = "*.meta.*.csv")
peppan_meta <- read.table(peppan_meta_read, header = TRUE, sep = "\t")
peppan_cluster_read <- list.files( path = "cluster_path", recursive = TRUE,
                                   full.names = TRUE, pattern = "*.cluster.csv")
peppan_cluster  <- read.table(peppan_cluster_read, header = TRUE, sep = "\t")
peppan_pagoo <- pagoo( data = roary_dataframe, org_meta = "roary_meta",
                       cluster_meta = "roary_cluster", sequences = 
                        "fasta_final", core_level = 100,
                       sep = "__")
panviz(peppan_presence_absence, location = "./roary_processing")
sequences <- roary_pagoo[["align"]]
write.FASTA(sequences, file = "./peppan_processing/selected_genes.fasta",
            header = TRUE)
alignment <- read.FASTA(file = "./peppan_processing/selected_genes.fasta",
                        type = "DNA")
alignment <- msa(alignment)
write.FASTA(alignment, file = "./peppan_processing/selected_genes_aligned.fasta",
            header = TRUE)
selected_phylogeny <- system2(command = "iqtree", args =
                               c("-s", "selected_genes_aligned.fasta" ,
                                 "-p", "selected_genes.phy" ,"-m" , "MFP+MERGE",
                                 "-B" , "1000", sep = " "))
outlier <- odseq(alignment, distance_metric = "affine", B = 1000,
                 threshold = 0.025)
alignment_plot <- ggmsa( alignment, start, end, color = "Shapely_NT",
                         font = "DroidSansMono", char_width = 0.5,
                         seq_name = TRUE) + geom_seqlogo(color = "Shapely_NT")
core_sequences <- roary_pagoo$core_sequences_4_phylogeny
(max_per_org = 1, fill = TRUE )
write.FASTA(core_sequences, file = "core_genes_peppan.fasta", header = TRUE)
read_core_genes <- read.FASTA(file = "core_genes_peppan.fasta", type = "DNA")
core_alignment <- msa(read_core_genes)
write.FASTA(core_alignment, file = "core_gene_peppan_aligned.fasta",
            header = TRUE)
core_phylogeny <- system2(command = "iqtree", args = c("-s",
                                     "core_gene_peppan_aligned.fasta" ,
                               "-p", "selected_genes.phy" ,"-m" , "MFP+MERGE",
                                          "-B" , "1000", sep = " "))
core_outlier <- odseq(core_alignment, distance_metric = "affine", B = 1000,
                      threshold = 0.025)
core_alignment_plot <- ggmsa( core_alignment, start, end, color = "Shapely_NT",
                              font = "DroidSansMono", char_width = 0.5,
                              seq_name = TRUE) + geom_seqlogo(color = "Shapely_NT")
+geom_msaBar() + geom_GC()
save.image(file = "peppan_alignment_plot", version = 1)
save.image(file = "core_roary_alignment_plot", version = 1)
pagoo_items <-  c( "pan_matrix", "organism", "clusters", "genes", 
                   "sequences","core_level", "core_genes", "core_clusters", 
                   "core_sequences", "cloud_genes","cloud_clusters", 
                   "cloud_sequences", "shell_genes", "shell_clusters", 
                   "summary_stats")
for (i in seq_along(pagoo_items))
{
 write.csv(peppan_pagoo$pagoo_items[i], file =
            "roary_pagoo$pagoo_items[i].csv", sep = ",")
}
peppan_pca <- peppan_pagoo$pan_pca()
peppan_pagoo$gg_pca(color = "organisms", size = 4) + theme_bw(base_size = 10) 
+scale_color_brewer("Set2") + geom_point() + 
distances = c("bray", "jaccard")
for (i in seq_along(distances)) {
 distances[i] <- roary_pagoo$dist(method = distances[i], binary = TRUE, 
                                  diag = TRUE)
}
bar_plot <- peppan_pagoo$gg_barplot()
bray_heatmap <- peppan_pagoo$ggdist(method = "bray")
jaccard_heatmap <- peppan_pagoo$ggdist(method = "jaccard")
peppan_pangenome_curves <- peppan_pagoo$gg_curves(what = "pangenome", size = 1) +
 ggtitle("pangenome curve of pirate_pagoo") + 
 geom_point(alpha = 0.1, size = 4) +
 theme_classic() + ylim(0,50000) + scale_color_brewer("Accent")
peppan_coregenome_curves <- peppan_pagoo$gg_curves(what = "coregenome", size = 1) +
 ggtitle("pangenome curve of pirate_pagoo") + 
 geom_point(alpha = 0.1, size = 4) +
 theme_classic() + ylim(0,50000) + scale_color_brewer(palette = "Set1") + 
 scale_color_gradient() + scale_fill_brewer(direction = -1)
peppan_pagoo$write_pangenome( dir = "peppan_processing", force = FALSE)
save.image(file = "bar_plt", version = 1)
save.image(file = "bray_heatmap", version = 1)
save.image(file = "jaccard_heatmap", version = 1)
save.image(file = "peppan_curves", version = 1)
#####################################################################
#################### PIRATE run #####################################
#####################################################################
if (argv$PIRATE_run) {
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
align <- argv$align
system("cp -r gff_path/*.gff /pirate_analysis")
if (length(list.files(path = "pirate_analysis/",
                      pattern = "*.gff" == length(list.files(path =
                                    "gff/path/", pattern = "*.gff")))))
{
 print("the gff files have been moved for the roary analysis
         and the file count has been verified")
}
system2(command = "PIRATE", args = c("-i", "./pirate_processing", "-a", 
                           "-r", "-k", "--diamond", "-s", "90,91,92,93,94,95", 
                           "--cd-step 1", "--cd-low 95"))
system2(command = "convert_to_roary.pl", args = 
                         "-i", "./pirate_processing/PIRATE.*.tsv", 
                                              "-o", "./pirate_to_roary.tsv")
pirate_dataframe <- read.table("./pirate_processing/pirate_to_roary.tsv", 
                                            header = TRUE, sep = "\t")
pirate_fasta <- list.files(path = "fasta_path", recursive = TRUE, 
                          full.names = TRUE, pattern = "*.fasta")
pirate_final <- sapply(roary_fasta, readDNAStringSet)
pirate_meta_read <- list.files( path = "meta_path", recursive = TRUE, 
                               full.names = TRUE, pattern = "*.meta.*.csv")
pirate_meta <- read.table(roary_meta_read, header = TRUE, sep = "\t")
pirate_cluster_read <- list.files( path = "cluster_path", recursive = TRUE,
                                  full.names = TRUE, pattern = "*.cluster.csv")
pirate_cluster  <- read.table(roary_cluster_read, header = TRUE, sep = "\t")
pirate_pagoo <- pagoo( data = roary_dataframe, org_meta = "roary_meta",
                      cluster_meta = "roary_cluster", sequences = "fasta_final", 
                      core_level = 100, sep = "__")
panviz(roary_presence_absence, location = "./roary_processing")
sequences <- pirate_pagoo[["align"]]
write.FASTA(sequences, file = "./pirate_processing/selected_genes.fasta",
            header = TRUE)
alignment <- read.FASTA(file = "./pirate_processing/selected_genes.fasta",
                        type = "DNA")
alignment <- msa(alignment)
write.FASTA(alignment, file = "./pirate_processing/selected_genes_aligned.fasta",
            header = TRUE)
selected_phylogeny <- system2(command = "iqtree", args =
                                    c("-s", "selected_genes_aligned.fasta" ,
                                 "-p", "selected_genes.phy" ,"-m" , "MFP+MERGE",
                                 "-B" , "1000", sep = " "))
outlier <- odseq(alignment, distance_metric = "affine", B = 1000,
                 threshold = 0.025)
alignment_plot <- ggmsa( alignment, start, end, color = "Shapely_NT",
                         font = "DroidSansMono", char_width = 0.5,
                         seq_name = TRUE) + geom_seqlogo(color = "Shapely_NT")
core_sequences <- roary_pagoo$core_sequences_4_phylogeny
(max_per_org = 1, fill = TRUE )
write.FASTA(core_sequences, file = "core_genes_pirate.fasta", header = TRUE)
read_core_genes <- read.FASTA(file = "core_genes_pirate.fasta", type = "DNA")
core_alignment <- msa(read_core_genes)
write.FASTA(core_alignment, file = "core_gene_pirate_aligned.fasta",
            header = TRUE)
core_phylogeny <- system2(command = "iqtree", args = c("-s",
                    "core_gene_pirate_aligned.fasta", "-p", "selected_genes.phy" 
                    ,"-m" , "MFP+MERGE", "-B" , "1000", sep = " "))
core_outlier <- odseq(core_alignment, distance_metric = "affine", B = 1000,
                      threshold = 0.025)
core_alignment_plot <- ggmsa( core_alignment, start, end, color = "Shapely_NT",
                              font = "DroidSansMono", char_width = 0.5,
                              seq_name = TRUE) + geom_seqlogo(color = "Shapely_NT")
+geom_msaBar() + geom_GC()
save.image(file = "roary_alignment_plot", version = 1)
save.image(file = "core_pirate_alignment_plot", version = 1)
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

pirate_pagoo$gg_pca(color = "organisms", size = 4) + theme_bw(base_size = 10) 
+scale_color_brewer("Set2") + geom_point() + facet_grid("~")
distances = c("bray", "jaccard")
for (i in seq_along(distances)) {
 distances[i] <- roary_pagoo$dist(method = distances[i], binary = TRUE, 
                                  diag = TRUE)
}
bar_plot <- pirate_pagoo$gg_barplot()
bray_heatmap <- pirate_pagoo$ggdist(method = "bray")
jaccard_heatmap <- pirate_pagoo$ggdist(method = "jaccard")
pirate_pangenome_curves <- pirate_pagoo$gg_curves(what = "pangenome", size = 1) +
 ggtitle("pangenome curve of pirate_pagoo") + 
 geom_point(alpha = 0.1, size = 4) +
 theme_classic() + ylim(0,50000) + scale_color_brewer("Accent")
pirate_coregenome_curves <- pirate_pagoo$gg_curves(what = "coregenome", size = 1) +
 ggtitle("pangenome curve of pirate_pagoo") + 
 geom_point(alpha = 0.1, size = 4) +
 theme_classic() + ylim(0,50000) + scale_color_brewer(palette = "Set1") + 
 scale_color_gradient() + scale_fill_brewer(direction = -1)
pirate_pagoo$write_pangenome( dir = "pirate_processing", force = FALSE)
save.image(file = "bar_plt", version = 1)
save.image(file = "bray_heatmap", version = 1)
save.image(file = "jaccard_heatmap", version = 1)
save.image(file = "pirate_curves", version = 1)
#####################################################################
############ Venn visualization of the common metagenomes ###########
#####################################################################
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
