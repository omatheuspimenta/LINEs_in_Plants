from Bio import SeqIO
from collections import defaultdict
from matplotlib.lines import Line2D
import matplotlib as plt
from pygenomeviz import GenomeViz
import pandas as pd

metadata_file = '../data/tandemrepeats_metadados_annotation.csv'
seqs_file = '../data/seq.fasta'
msat_file = '../results/csv/'
metadata = pd.read_csv(metadata_file, sep=';', index_col=False)
elements_files = metadata['Name of elements'].to_list()


def make_plot(seq_name:str,
              seq_len:int,
              minisats_regions: list[tuple[int,int]], 
              regions: defaultdict[list[tuple[int,int]]]):
    """
    Generates a plot of genomic regions with specified features and saves it as an image file.
    Args:
        seq_name (str): The name of the sequence to be plotted.
        seq_len (int): The length of the sequence.
        minisats_regions (list[tuple[int, int]]): A list of tuples representing the start and end positions of minisatellite regions.
        regions (defaultdict[list[tuple[int, int]]]): A defaultdict containing lists of tuples representing the start and end positions of various genomic regions, keyed by region name.
    Returns:
        None
    """
    
    colors = {'MicroSats': 'grey',
              'ORF1 RRM': 'cyan',
              'ORF1 DUF': 'red',
              'ORF1 Zf-CCHC': 'purple',
              'ORF1 PHA': 'skyblue',
              'ORF2 EEP': 'green',
              'ORF2 RT': 'blue',
              'ORF2 ZF-RVT': 'violet',
              'ORF2 RH': 'orange',
              '(GTT)n': 'crimson',
              'Poli-A': 'black'}
    handles = []
    
    gv = GenomeViz(fig_track_height=0.5, feature_track_ratio=2)
    # gv.set_scale_bar(ymargin=1.5)
    gv.set_scale_xticks()
    
    track = gv.add_feature_track(seq_name, seq_len)
    
    track.add_sublabel()
    
    # Add styled features
    for msat in minisats_regions:
        start, end = msat[0], msat[1]
        track.add_feature(start,
                          end,
                          -1,
                          plotstyle="box",
                          fc="grey")
                          # label="MicroSats",
                          # )#,
                          #text_kws=dict(rotation=-45, vpos="bottom"))
    handles.append(Line2D([], [],
                          marker="s",
                          color="grey",
                          label="MicroSats",
                          ms=12,
                          ls="none"))
    
    for k, v in regions.items():
        track.add_feature(v[0][0],
                          v[0][1],
                          1,
                          plotstyle="box",
                          fc=colors.get(k))
        handles.append(Line2D([], [],
                              marker="s",
                              color=colors.get(k),
                              label=k,
                              ms=12,
                              ls="none"))
    
    fig = gv.plotfig()
    
    # Plot legend for groups
    _ = fig.legend(
        handles=handles,
        fontsize=12,
        title="Groups",
        title_fontsize=12,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        handlelength=1.0,
    )
    
    seq_name = seq_name.split('.')[0]
    plot_name = '../TandemRepeats/results/plots/' + seq_name + 'png'
    fig.savefig(plot_name, dpi=fig.dpi)
    fig.clf()
    
def get_minisat_regions(df) -> list[tuple[int,int]]:
    """
    Extracts regions of minisatellites from a DataFrame.
    Args:
        df (pandas.DataFrame): A DataFrame containing 'start' and 'end' columns 
                                with integer values representing the start and end 
                                positions of minisatellite regions.
    Returns:
        list[tuple[int, int]]: A list of tuples where each tuple contains two integers 
                                representing the start and end positions of a minisatellite region.
    """
    
    
    minisats_regions = []
    
    for srt, end in zip(df['start'],df['end']):
        minisats_regions.append((int(srt),int(end)))

    return minisats_regions

def get_regions(df) -> defaultdict[list[tuple[int,int]]]:
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
            regions['ORF1 RRM'].append((int(row['ORF 1 RRM início']), int(row['ORF 1 RRM fim'])))
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
        if seq_name in elements_files:
            _msat_file = msat_file + seq_name + '.csv'
            temp_df = metadata[metadata['Name of elements'] == seq_name]
            try:
               msat_df = pd.read_csv(_msat_file, sep=';', index_col=False)
               print("running...", str(record.id))
               # seq_len = len(record)
               seq_len = int(temp_df['Size (bp)'].iloc[0])
               if seq_len < len(record):
                   seq_len = len(record)
               
               minisats_regions = get_minisat_regions(msat_df)
               regions = get_regions(temp_df)
               
               make_plot(seq_name, seq_len, minisats_regions, regions)
            except Exception as e:
                print("NO MICROSAT FOR", seq_name)
                print(e)
