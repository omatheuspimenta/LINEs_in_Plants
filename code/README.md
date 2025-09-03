# Using the `get_domains` validation script

To run the script, we strongly recommend using **`uv`**.  
If you don’t have `uv` installed, you can find installation instructions here: [https://docs.astral.sh/uv/](https://docs.astral.sh/uv/)

```bash
Usage: get_domains.py [-h] [--k K] [--s S] --input INPUT --metadata METADATA

Validate domain sequences using the k-mers approach.

Options:
  -h, --help           Show this help message and exit
  --k K                Length of k-mers
  --s S                Step size for k-mers
  --input INPUT        Path to the input directory containing domain files
  --metadata METADATA  Path to the metadata file
```

The `get_domains.py` script requires two mandatory parameters:

* **`--input`**: the folder containing domain sequences in FASTA format (one domain per FASTA file). Example:

```
├── EEP.fasta
├── ORF1 DUF.fasta
├── ORF1 Zf-CCHC.fasta
├── Poli-A.fasta
├── RH.fasta
├── RT.fasta
└── ZF-RVT.fasta
```

* **`--metadata`**: the metadata file, where each column corresponds to a domain in the following format:

```csv
Organism Name;SeqName;Name of element;Size (bp);ORF-1 RRM;ORF-1 DUF;ORF-1 zf-CCHC;EEP;RT;RVT;RH;GTT;Poli-A
Lunularia cruciata (L.) Dumort. ex Lindb.;RTE_Lunularia_cruciata_seq1;RTE_seq1;4313;;;;537-1328;2058-2822;;;;4305-4313
Lunularia cruciata (L.) Dumort. ex Lindb.;RTE_Lunularia_cruciata_seq2;RTE_seq2;3685;;;;512-1303;2033-2797;;;;3675-3685
```

To execute the script with `uv`, simply run the following in your terminal:

```bash
uv run get_domains.py --input input_folder --metadata metadata.csv
```