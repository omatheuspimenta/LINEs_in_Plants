# LINEs in Plants
## Table of Contents
1. [Introduction](#introduction)
2. [Repository Contents](#repository)
3. [Citation](#cite)

## Introduction

Long interspersed nuclear elements (LINEs) are a type of non-LTR retrotransposon belonging to the Retroposineae suborder. They are found in all eukaryotes. Compared to LTR retrotransposons, these elements are underrepresented in plant genomes and are also poorly understood in the wider context of land plants. In this study, we reorganised an RT-based database and developed strategies for searching for and filtering LINEs in 54 plant genomes. Analysis was performed on species from the classes Marchantiopsida, Sphagnopsida, Bryopsida, Polypodiopsida, Pinopsida and Magnoliopsida. We investigated the structure, evolution and abundance of LINEs in these genomes based on the conservation of probable autonomous elements. Two families of elements have dominated plant genomes: L1 and RTE. However, fragments of degenerate I, R2 and Jockey elements were also identified, although in low percentages. Most elements exhibited significant variability in ORF1 (5'), while ORF2 remained conserved within two 'essential' domains (EEP and RT), as well as two variable and non-essential domains (Zf-RVT and RNH). Only those containing a 3' poly-A tail were considered for the annotation of probable autonomous elements. Phylogenetic analysis of the reverse transcriptase sequence revealed two consistent clades: one containing elements of the RTE superfamily and one containing elements of L1. Each of these clades contained subclades.

## Repository
This repository is organized as follows:

### Code
The `code` folder contains all files and instructions related to the replication of the methodology:

### Data
The `data` folder contains all files related to the work and is subdivided as follows:

- [`metadata.csv`](data/metadata.csv): Contains the metadata used for annotation, extraction, and database generation.
- [`metadata.xlsx`](data/metadata.xlsx): Contains the metadata used for annotation, extraction, and database generation.
- [`annotation.xlsx`](data/annotation.xlsx): Contains the metadata used for annotation, extraction, and database generation.
- [`LINE_domain_annotations.gff3`](data/LINE_domain_annotations.gff3): Conserved domains annotated by elements and species (GFF3)
- [`All_LINEs_EEP_RVT_RH_prot_DB_RepBase_REXdb_NCBI.fasta`](data/All_LINEs_EEP_RVT_RH_prot_DB_RepBase_REXdb_NCBI.fasta): Contains sequences of LINE elements Exonuclease-Endonuclease-Phosphatase (EEP), zinc-binding region in Reverse Transcriptase (zf-RVT) and Ribonuclease-H (RH) domains that were present in online databases (RepBase, REXdb and NCBI).  
- [`All_LINES_nucl_complete_sequences.fa`](data/All_LINES_nucl_complete_sequences.fa): Database containing the complete nucleotide sequences of all annotated LINEs elements.  
- [`All_LINEs_RT_prot_DB_RepBase_REXdb_NCBI.fasta`](data/All_LINEs_RT_prot_DB_RepBase_REXdb_NCBI.fasta): Contains sequences of LINE elements Reverse Transcriptase (RT) domains that were present in online databases (RepBase, REXdb and NCBI).  
- [LINEs_EEP_DB_nucl_extracted_seqs.fasta](data/LINEs_EEP_DB_nucl_extracted_seqs.fasta): Database containing nucleotide sequences of the Exonuclease-Endonuclease-Phosphatase (EEP) domains found in the annotated LINEs elements.  
- [LINEs_EEP_DB_prot_extracted_seqs.fasta](data/LINEs_EEP_DB_prot_extracted_seqs.fasta): Database containing amino acid sequences of the Exonuclease-Endonuclease-Phosphatase (EEP) domains found in the annotated LINEs elements.  
- [LINEs_RH_DB_nucl_extracted_seqs.fasta](data/LINEs_RH_DB_nucl_extracted_seqs.fasta): Database containing nucleotide sequences of the Ribonuclease-H (RH) domains found in the annotated LINEs elements.  
- [LINEs_RH_DB_prot_extracted_seqs.fasta](data/LINEs_RH_DB_prot_extracted_seqs.fasta): Database containing amino acid sequences of the Ribonuclease-H (RH) domains found in the annotated LINEs elements.  
- [LINEs_RT_DB_nucl_extracted_seqs.fasta](data/LINEs_RT_DB_nucl_extracted_seqs.fasta): Database containing nucleotide sequences of the Reverse Transcriptase (RT) domains found in the annotated LINEs elements.  
- [LINEs_RT_DB_prot_extracted_seqs.fasta](data/LINEs_RT_DB_prot_extracted_seqs.fasta): Database containing amino acid sequences of the Reverse Transcriptase (RT) domains found in the annotated LINEs elements.  
- [LINEs_RVT_DB_nucl_extracted_seqs.fasta](data/LINEs_RVT_DB_nucl_extracted_seqs.fasta): Database containing nucleotide sequences of the zinc-binding region in Reverse Transcriptase (zf-RVT) domains found in the annotated LINEs elements.  
- [LINEs_RVT_DB_prot_extracted_seqs.fasta](data/LINEs_RVT_DB_prot_extracted_seqs.fasta): Database containing amino acid sequences of the zinc-binding region in Reverse Transcriptase (zf-RVT) domains found in the annotated LINEs elements.  


Additionally, the `data` folder contains the following folder: `domains`. This folder contain all sequences, domains and figures. Each organism has its own folder with the following structure:

- `Annotation`: Contains annotations obtained from CD-HIT and microsat annotations.
- `LINES`: Contains the complete LINE sequences for the organism.
- `figs`: Contains the LINE figures obtained through the annotation and machine learning validation approach.

## Cite
Soon
