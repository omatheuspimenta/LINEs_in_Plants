# LINEs in Plants
## Table of Contents
1. [Introduction](#introduction)
2. [Repository Contents](#repository)
3. [Citation](#cite)

## Introduction

Long interspersed nuclear elements (LINEs) are an order of non-LTR retrotransposons found in all eukaryotes. These transposable elements are more abundant in animals than in plants and are poorly understood in terms of structure, ORF conservation, classification and evolution. In this study, we reorganized a reverse transcriptase (RT)-based database and developed a search and filtering strategy for LINEs in 40 plant genomes sequenced by long molecule technology. We performed a comparative and systematic analysis of LINEs in species of Hepatophyta, Bryophyta, Polypodiophyta, Pinophyta, and Angiospermae (eudicots and monocots). Aspects of structure, evolution, dating, abundance, and distribution of LINEs along chromosomes were evaluated.Based on the conservation of intact LINEs, we found that three lineages predominate in plants: I, L1, and RTE.  In most elements, ORF1 was located at the 5' end, but they were highly variable in the presence of the three subdomains (RRM, DUF and zf-CHCC) because they were truncated or absent. In contrast, ORF2 was composed of two conserved and essential domains, exo-endonuclease phosphatase (EEP) and reverse transcriptase (RT), among other non-essential domains, (Zf-RVT and RNH). All annotated LINEs presented one or more A-rich sequences in the terminal 3' region (poly-A tail), but the flanking UTR regions were not successfully identified due to the immense variation and difficulty of annotation. Comparative analysis of LINE sequences showed that the set of intact and degenerate elements represented up to 3% of the genomes, with no clear relationship between the accumulation of LINEs and genome size. I elements were the largest, with the most gene domains and expansion peaks at ~10 to 15 million years ago. They were followed by L1 and RTE elements. The I superfamily predominated in the most basal plant groups, whereas RTE elements (present in monocots) were the smallest, with only EEP and RT domains and expansion peaks at ~5 my. The overall size of LINEs was quite variable, especially in the regions between gene domains, which accounted for an average of 58% of the total size of the elements. The distribution patterns of LINEs in the pseudochromosomes of plant genomes showed great heterogeneity in distribution and accumulation along chromosomes, with fragmented elements being more abundant than whole ones. The alternative annotation of LINEs using a machine learning tool based on short K-mers efficiently recognized all gene domains and showed that these interstitial regions have signatures that can also be useful for detecting and identifying LINEs quickly and effectively in systematic analyses. The analyses suggest that LINEs in plants have undergone evolutionary processes involving breakage and loss of regions and domains, such that I elements would be the most basal in the plant LINE phylogeny and RTE elements the most derived. In addition to being fast and using few hardware resources, ML was more efficient than traditional alignment in detecting small and large changes involving essential and non-essential domains, and helped us understand the degree of conservation of complete elements.

## Repository
This repository is organized as follows:

### Code
The `code` folder contains all files related to the replication of the methodology:

- [`compress.sh`](code/compress.sh): Compresses all FASTA files to `.tar.xz` format.
- [`get_kmers.ipynb`](code/get_kmers.ipynb): Generates the RAW, soft, and hard annotation regions, along with annotation images.
- [`organize_files_domains.py`](code/organize_files_domains.py): Organizes the LINE `.fasta` files into `.fasta` files of the LINE domains.
- [`organize_files_lines.py`](code/organize_files_lines.py): Organizes the LINE `.fasta` files.
- [`organize_files_organisms.py`](code/organize_files_organisms.py): Organizes the LINE `.fasta` files of the organisms into separate folders.
- [`plot_regions_msat.py`](code/plot_regions_msat.py): Generates graphs of the annotations, including microsat annotations.
- [`run_RPTRF.py`](code/run_RPTRF.py): Runs the RPTRF program and organizes and analyzes the results.
- [`split_seqs_regions.py`](code/split_seqs_regions.py): Separates the regions of the LINEs domains based on the metadata.

### Data
The `data` folder contains all files related to the work and is subdivided as follows:

- [`Article_Figures_Tables`](data/Article_Figures_Tables/): Contains all figures and tables related to the article, further divided into:
  - [`Article_Figures`](data/Article_Figures_Tables/Article_Figures/): Contains all figures present in the article.
  - [`Supplementary_Figures`](data/Article_Figures_Tables/Supplementary_Figures/): Contains all supplementary figures of the article.
  - [`Supplementary_Tables`](data/Article_Figures_Tables/Supplementary_Tables/): Contains all supplementary tables of the article.
- [`Database_RT_LINEs`](data/Database_RT_LINEs/): Contains the UGENE Kimura parameters, the initial database, and all plant RTs.
- [`anotacao_graph.csv`](data/anotacao_graph.csv): Contains the annotations used to generate supplementary images.
- [`metadata_updated.csv`](data/metadata_updated.csv): Contains the metadata used for annotation, extraction, and database generation.
- [`microsats_analysis.csv`](data/microsats_analysis.csv): Contains the microsat analyses of all sequences after the execution of the RPTRF program.
- [`seqs.tar.xz`](data/seqs.tar.xz): Contains the sequences of the complete LINEs used for the analyses.
- [`tandemrepeats_metadados_annotation.csv`](data/tandemrepeats_metadados_annotation.csv): Contains the metadata used for microsat analyses.

Additionally, the `data` folder contains the following subfolders: `Angiospermophyta`, `Bryophyta`, `Hepatophytas`, `Pinophyta`, and `Polypodiophyta`. These folders contain all sequences, annotations, domains, figures, and analyses. Each organism has its own folder with the following structure:

- `Annotation`: Contains annotations obtained from CD-HIT and microsat annotations.
- `Entire_LINES_sequences`: Contains the complete LINE sequences for the organism.
- `domains`: Contains domain sequences separated into `Separate_sequence_Domains_LINEs` and combined sequences in `Combined_sequence_Domains_LINEs`.
- `Test_Sequences`: When available, this folder contains sequences without annotations that were used to validate the proposed k-mer methodology.

## Cite
Soon