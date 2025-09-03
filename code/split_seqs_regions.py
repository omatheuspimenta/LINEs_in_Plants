from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
from collections import defaultdict
from pathlib import Path

metadata_file = r'../data/metadata.csv'
seqs_file = r'../data/All_LINES_nucl_complete_sequences.fa'
output_folder = r'../data/domains/'
metadata = pd.read_csv(metadata_file, sep=';', index_col=False)
elements_files = metadata['SeqName'].to_list()

def get_regions(df) -> defaultdict(list):
    """
    Extracts specific regions from a DataFrame and organizes them into a defaultdict of lists.
    This function iterates over each row of the provided DataFrame and checks for non-null values
    in specific columns that represent the start and end positions of various regions. It then 
    appends these positions as tuples to the corresponding lists in the defaultdict.
    Args:
        df (pandas.DataFrame): The input DataFrame containing columns with region start and end positions.
    Returns:
        defaultdict(list): A dictionary where keys are region names and values are lists of tuples,
                           each tuple containing the start and end positions of the region.
    """
    
    # mapping of df column → regions key
    col_map = {
        "ORF-1 RRM": "ORF1 RRM",
        "ORF-1 DUF": "ORF1 DUF",
        "ORF-1 zf-CCHC": "ORF1 Zf-CCHC",
        "EEP": "EEP",
        "RT": "RT",
        "RVT": "ZF-RVT",
        "RH": "RH",
        "GTT": "(GTT)n",
        "Poli-A": "Poli-A"
    }
    
    regions = defaultdict(list)
    
    for col, key in col_map.items():
        # take only non-null values
        valid = df[col].dropna()
    
        if valid.empty:
            continue  # nothing to do for this column
    
        # split into start / end with regex (safer than plain split)
        split_vals = (
            valid.str.extract(r"(\d+)-(\d+)")
                 .dropna()
                 .astype(int)
        )
    
        regions[key].extend(list(zip(split_vals[0], split_vals[1])))
        
    return regions



with open(seqs_file, encoding="utf-8") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        seq_name = str(record.id)
        seq = str(record.seq).upper()
        if seq_name in elements_files:
            temp_df = metadata[metadata['SeqName'] == seq_name]
            
            # if (len(seq) + 1) != temp_df['Size (bp)'].iloc[0]:
            #     print(seq_name, ": Seq_Len:", len(seq) + 1, "Seq_Len (annotated):",  temp_df['Size (bp)'].iloc[0])
        # else:
            # print(seq_name)
            regions = get_regions(temp_df)
            for k, v in regions.items():
                start = v[0][0] - 1
                stop = v[0][1]
                cutted = ''.join(list(seq)[start:stop])
                # if len(cutted) == 0:
                #     print(seq_name, " ", k)
                seq_file_name = seq_name + '_' + str(k) + '_' + str(start+1) + ':' + str(stop)
                cutted_seq = SeqRecord(
                    Seq(cutted),
                    id=seq_file_name,
                    description=''
                )
                
                seq_ = SeqRecord(
                    Seq(seq),
                    id=seq_name,
                    description=''
                )
                
                organism = temp_df['Organism Name'].iloc[0]  # get the first value
                save_path = Path(output_folder) / organism
                save_path_lines = save_path / 'LINES'
                
                save_path.mkdir(parents=True, exist_ok=True)
                save_path_lines.mkdir(parents=True, exist_ok=True)
                                
                file_name = str(k) + '.fasta'
                output_file = save_path / file_name
                line_name = seq_name + '_LINE.fasta'
                output_file_line = save_path_lines / line_name
                
                with open(file=output_file, mode = "a") as out_file:
                    SeqIO.write(cutted_seq, out_file, "fasta")
                    
            with open(file=output_file_line, mode = "a") as out_file:
                SeqIO.write(seq_, out_file, "fasta")