# MetaPanRan
MetaPanRan

An integrated pipeline for the prediction of the bacterial genomes and analysis and visualization of the pangenomes from ROARY, PIRATE, PEPPAN, PANROO. This allows you to run the three accurate pangenome predictions from ROARY, PEPPAN, PIRATE, PANROO and then gives the pangenomes. It also gives a comparative assessment of the pangenomes
that is conserved across all through a venn visualization. There are several versions and is in beta testing and the complete package will be updated with the final release as MetaPanRan_release.R

However, you are free to try it and report the bugs and i would be more than happy to acknowledge you in the git contributions. 

A brief manuscript version is present below for reading as a manual. References to be updated. 

MetaPanRan: A standalone R based streamline workflow for analyzing metagenomes from pangenomics to machine learning. 
Gaurav Sablok, Thulani Peter Makhalanyane
Department of Biochemistry, Genetics and Microbiology,
Faculty of Natural and Agricultural Sciences, 
UNIVERSITY OF PRETORIA

Abstract: 
 Owing to the wide diversity and growing enormous potential of metagenomics widely associated with many of the biological spans, it is detrimental to have potential look at the meta pangenomics. In the present paper, we present an integrated work-flow for analyzing the microbial pangenomes using several state of the art pangenomics methods such as ROARY, PANAROO, PIRATE and PEPPAN. MetaPanRan allows you to run all these methods as a part of the application and allows for the downstream analysis of the microbial pangenomes such as the alignments, core alignments, alignment filtering, phylogeny, estimation of dn/ds, detection of outliers, and visualization and estimation of pangenomes. MetaPanRan serve as a common platform for the detection of the metagenomics and downstream analysis associated with metagenomes. Microbes play an important role in shaping the diversity and classification of the organisms as a result reflecting changes in the health as well as the diversity. Application of next generation sequencing methods has not only widely elaborated the state of the art metagenomics by illustrating individual microbial strains as well as the hierarchy level of the metagenomics at the strain level. Recently, a surge of metagenomics application has seen in understanding the pangenomics of the microbes, which has allowed for identifying common traits of interest as well as parallel application of genomic breadth of approaches at a multi-level class and at the phylum level. Several state of the art tools have been developed for the identification of the metagenomes such as ROARY (Reference), PANAROO (Reference), PIRATE (Reference), PEPPAN (Reference) and individual packages such as micropan (Reference), Pagoo (Reference), Panstripe (Reference) have been developed for the downstream analysis. However, the common limitation stays with the non-availability of the streamline work-flow for analyzing the pangenomes all togetherIn this paper, we present MetaPanRan which addresses this problem by providing a streamline work-flow by integrating all these metagenomes prediction tools. In addition, MetaPanRan also provide options for doing the genome predictions using the Prokka(Reference) and can be combined with the additional provided user supplied annotation files (format = “*.gff”) for the downstream pangenome predictions and annotations. MetaPanRan work-flow is depicted in Figure 1. It allows the users to predict the meta-pangenomes using the ROARY (Reference), PANAROO (Reference), PIRATE (Reference), PEPPAN (Reference). Following the prediction of the metagenomes, it performs several downstream analysis such as doing the pangenome alignments, alignment visualization, masking the alignment for the gaps, and implementing a model based phylogeny on the selected core genes and the pangenome alignments. It uses muscle for the alignment and also for the staggered alignment and IQTREE(Reference) under the “MFP+MERGE” model with 1000 bootstraps for the clade phylogeny. Maximum likelihood phylogeny has been the main stream of development of the phylogeny with or without constraints. MetaPanRan implements the maximum likelihood based phylogeny using TreeLine with a cutoff of 0.5 (Reference)Sequence visualization and looking at the sequence plots with the base probabilities is one of the main stream of looking at the consensus alignments. MetaPanRan implements ggmsa (Reference) for the selection and visualization of the base probabilities. Panstripe (Reference) has been integrated in MetaPanRan for the visualization of the core pan phylogeny and plotting the “gaussian” distribution of the pan metagenomes. MetaPanRan requires minimal installation and  it performs all the installation and directory wise hierarchical classification as a part of the pan genome predictions and downstream analysis. All the core and pan genes files as a part of the user selected genes are written as fasta formatted files and base probabilities. 
 Microbial genes shows accelerated rate of evolution and this has been widely demonstrated and has been used for simulation and genetic algorithms development. To understand the accelerated rate of evolution, MetaPanRan first estimates the alignment using Muscle(Reference) and then delete cols with gaps using the APE package(Reference) and subsequently detected the outliers in the alignment using OdSeq(Reference). Following the outlier detection and alignment filtering, it estimates the transition/transversion rates using the seqinr package(Reference) and write a dataframe object. To conclude, we present an integrated version of metagenome predictions, analysis and visualization as a stand alone R tool, which allows the users to do meta pangenome analysis with minimal interactions and at a fine scale of the metagenome. This bridges a lack of the integrated workflow for the meta pangenomes. 


References:
1.	Protocol for post-processing of bacterial pangenome data using Pagoo pipeline. STAR Protocols 2, 100802, December 17, 2021
2.	Snipen and Liland BMC Bioinformatics (2015) 16:79 micropan: an R-package for microbial pan-genomics 
3.	microshades: An R Package for Improving Color Accessibility and Organization of Microbiome Data
4.	Schliep, K.P. (2011). phangorn: phylogenetic analysis in R. Bioinformatics 27, 592–593.
5.	Seemann, T. (2014). Prokka: rapid prokaryotic genome annotation. Bioinformatics 30, 2068– 2069.
6.	R.A., et al. (2020). Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biol. 21, 180.
7.	Wright, E.S. (2015). DECIPHER: harnessing local sequence context to improve protein multiple sequence alignment. BMC Bioinformatics 16, 322.
8.	Zhou, Z., Charlesworth, J., and Achtman, M. (2020). Accurate reconstruction of bacterial pan- and core genomes with PEPPAN. Genome Res. 30, 1667– 1679.
9.	Bayliss, S.C., Thorpe, H.A., Coyle, N.M., Sheppard, S.K., and Feil, E.J. (2019). PIRATE: A fast and scalable pangenomics toolbox for clustering diverged orthologues in bacteria. GigaScience 8, giz119.
10.	MicrobiotaProcess: A comprehensive R package for deep mining microbiome

Gaurav Sablok \
Senior Postdoctoral Fellow \
Faculty of Natural and Agricultural Sciences Room 7-35, \
Agricultural Sciences Building University of Pretoria, \
Private Bag X20 Hatfield 0028, South Africa \
ORCID: https://orcid.org/0000-0002-4157-9405 \
WOS: https://www.webofscience.com/wos/author/record/C-5940-2014 \
Github:https://github.com/sablokgaurav \
Linkedln: https://www.linkedin.com/in/sablokgaurav/ \
ResearchGate: https://www.researchgate.net/profile/Gaurav-Sablok \
Academia : https://up-za.academia.edu/GauravSablok \
Plant Bioinformatics | Illumina | PacBio | OxfordNanoPore | Transcriptomics | \
Genomics | Metagenomics | DevOPS | Data Engineering and Analyst | \
Python | R | Java | Scala | Bash | DJANGO | Machine Learning \

