from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
from collections import defaultdict

metadata_file = r'../data/metadata_updated.csv'
seqs_file = r'../data/seqs.fasta'
output_folder = r'../data/domains/'
metadata = pd.read_csv(metadata_file, sep=';', index_col=False)
elements_files = metadata['Name of elements'].to_list()

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
    
    regions = defaultdict(list)

    # Iterate over rows to access individual values
    for idx, row in df.iterrows():
        if pd.notna(row['ORF 1 RRM início']):
            regions['ORF1 RRM'].append((int(row['ORF 1 RRM início']), int(row['ORF 1 RRM  fim'])))
        if pd.notna(row['ORF 1 DUF início']):
            regions['ORF1 DUF'].append((int(row['ORF 1 DUF início']), int(row['ORF 1 DUF fim'])))
        if pd.notna(row['ORF 1 Zf-CCHC início']):
            regions['ORF1 Zf-CCHC'].append((int(row['ORF 1 Zf-CCHC início']), int(row['ORF 1 Zf-CCHC fim'])))
        if pd.notna(row['ORF 1 PHA(?) início']):
            regions['ORF1 PHA'].append((int(row['ORF 1 PHA(?) início']), int(row['ORF 1 PHA(?)  fim'])))
        if pd.notna(row['ORF 2 EEP início']):
            regions['ORF2 EEP'].append((int(row['ORF 2 EEP início']), int(row['ORF 2 EEP  fim'])))
        if pd.notna(row['ORF 2 RT início']):
            regions['ORF2 RT'].append((int(row['ORF 2 RT início']), int(row['ORF 2 RT fim'])))
        if pd.notna(row['ORF 2 Zf-RVT início']):
            regions['ORF2 ZF-RVT'].append((int(row['ORF 2 Zf-RVT início']), int(row['ORF 2 Zf-RVT  fim'])))
        if pd.notna(row['ORF 2 RH início']):
            regions['ORF2 RH'].append((int(row['ORF 2 RH início']), int(row['ORF 2 RH fim'])))
        if pd.notna(row['(GTT)n início']):
            regions['(GTT)n'].append((int(row['(GTT)n início']), int(row['(GTT)n fim'])))
        if pd.notna(row['Poli-A início']):
            regions['Poli-A'].append((int(row['Poli-A início']), int(row['Poli-A fim'])))
    
    return regions



with open(seqs_file, encoding="utf-8") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        seq_name = str(record.id)
        seq = str(record.seq).upper()
        if seq_name in elements_files:
            temp_df = metadata[metadata['Name of elements'] == seq_name]
        # else:
            # print(seq_name)
            regions = get_regions(temp_df)
            for k, v in regions.items():
                start = v[0][0] - 1
                stop = v[0][1]
                cutted = ''.join(list(seq)[start:stop])
                
                seq_file_name = seq_name + '_' + str(k) + '_' + str(start+1) + ':' + str(stop)
                cutted_seq = SeqRecord(
                    Seq(cutted),
                    id=seq_file_name,
                    description=''
                )
                
                output_file = output_folder + str(k) + seq_name + '.fasta'
                
                with open(file=output_file, mode = "w") as out_file:
                    SeqIO.write(cutted_seq, out_file, "fasta")
            