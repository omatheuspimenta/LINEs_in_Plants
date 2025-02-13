from collections import Counter, defaultdict
from itertools import islice
from Bio import SeqIO
import subprocess
import pandas as pd

def get_sequence_interval(row: pd.Series) -> pd.Interval:
    """
    Create a closed interval using the 'start' and 'end' values from a row of data.

    This function constructs a closed interval using the 'start' and 'end' values from
    a pandas Series row. The interval is closed on both ends.

    Args:
        row (pd.Series): A pandas Series containing 'start' and 'end' values.

    Returns:
        pd.Interval: A closed interval defined by the 'start' and 'end' values.
    """
    return pd.Interval(row.start, row.end, closed='both')


file_in = 'data/seqs.fasta'
file_temp = 'temp.fasta'
file_result = 'result-1.txt'
output_file = 'analysis.csv'
command = './RPTRF -s temp.fasta -m 20 -t 5'
col_names = ['start', 'end', 'len', 'motif_seq', 'seq']
col_names_output = ['seq_name', 'seq_len', 'n_microsat', 'microsat_sum', 'microsat_ratio',
                    'microsat_position_begin', 'microsat_position_end',
                    'microsat_A', 'microsat_C', 'microsat_G', 'microsat_T',
                    'microsat_AC', 'microsat_AG', 'microsat_AT',
                    'microsat_CG', 'microsat_CT',
                    'microsat_GT']

with open(output_file, mode='a', encoding='utf-8') as out_file:
    out_file.write(';'.join(col_names_output))
    out_file.write('\n')

with open(file_in, encoding="utf-8") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        SeqIO.write(record, file_temp, 'fasta')
        print("running...", str(record.id))
        subprocess.run(command,
                       shell=True,
                       check=True,
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.STDOUT)
        
        subprocess.run('rm temp.fasta',
                       shell= True, 
                       check=True)
        
        seq_len = len(record)
        seq_middle = seq_len//2
        
        first_region = pd.Interval(left=0, right=seq_middle, closed='both')
        last_region = pd.Interval(left=seq_middle, right=seq_len, closed='both')
        
        values = []
        values_output = []
        with open(file_result, encoding="utf-8") as result_file:
            for line in islice(result_file, 5, None):
                values.append(line.split())
        
        # subprocess.run('rm ' + file_result, 
        #                shell=True, 
        #                check=True)
        subprocess.run('mv ' + file_result + ' results/txt/' + str(record.id) + '.txt', 
                       shell=True, 
                       check=True)
        
        if len(values) > 0:
            
            df = pd.DataFrame(values, columns=col_names)
            df['start'] = df['start'].astype(int)
            df['end'] = df['end'].astype(int)
            df['len'] = df['len'].astype(int)
            df['motif_seq'] = df['motif_seq'].astype(str)
            df['seq'] = df['seq'].astype(str)
            microsat_interval = df.apply(get_sequence_interval, axis=1)
            
            microsat_sum = df['len'].sum()
            
            values_output.append(record.id)
            values_output.append(str(seq_len))
            values_output.append(str(df.shape[0]))
            values_output.append(str(microsat_sum))
            values_output.append(str(microsat_sum/seq_len))
            
            first_region_count = 0
            last_region_count = 0
            for interval in microsat_interval:
                if interval in first_region:
                    first_region_count += 1
                elif interval in last_region:
                    last_region_count += 1
            
            values_output.append(str(first_region_count))
            values_output.append(str(last_region_count))
            
            microsat_nuc = defaultdict(int)
            for element in df['seq']:
                counter_temp = Counter(element)
                if len(counter_temp) > 1 and counter_temp[max(counter_temp)] - counter_temp[min(counter_temp)] > 0:    
                    microsat_nuc[max(counter_temp)] += 1
                if len(counter_temp) == 1:
                    microsat_nuc[max(counter_temp)] += 1
            
            values_output.append(str(microsat_nuc['A']))
            values_output.append(str(microsat_nuc['C']))
            values_output.append(str(microsat_nuc['G']))
            values_output.append(str(microsat_nuc['T']))
            
            at, ag, ac, tg, tc, gc = 0, 0, 0, 0, 0, 0
            for element in df['motif_seq']:
                temp = element.split('(')
                if temp[0] == '2':
                    temp = temp[1].split(')')
                    if (temp[0][0]=="A" and temp[0][1]=="T") or (temp[0][0]=="T" and temp[0][1]=="A"):
                        at += 1
                    if (temp[0][0]=="A" and temp[0][1]=="G") or (temp[0][0]=="G" and temp[0][1]=="A"):
                        ag += 1
                    if (temp[0][0]=="A" and temp[0][1]=="C") or (temp[0][0]=="C" and temp[0][1]=="A"):
                        ac += 1
                    if (temp[0][0]=="T" and temp[0][1]=="G") or (temp[0][0]=="G" and temp[0][1]=="T"):
                        tg += 1
                    if (temp[0][0]=="T" and temp[0][1]=="C") or (temp[0][0]=="C" and temp[0][1]=="T"):
                        tc += 1
                    if (temp[0][0]=="G" and temp[0][1]=="C") or (temp[0][0]=="C" and temp[0][1]=="G"):
                        gc += 1
                
            values_output.append(str(ac))
            values_output.append(str(ag))
            values_output.append(str(at))
            values_output.append(str(gc))
            values_output.append(str(tc))
            values_output.append(str(tg))
            
            with open(output_file, mode='a', encoding='utf-8') as out_file:
                out_file.write(';'.join(values_output))
                out_file.write('\n')
        else:
            values_output = ['NO_MICROSAT_FOUND'] * 17
            values_output[0] = str(record.id)
            with open(output_file, mode='a', encoding='utf-8') as out_file:
                out_file.write(';'.join(values_output))
                out_file.write('\n')
        
        df.to_csv('results/csv/' + str(record.id) + '.csv', sep=';', index=False)
        print(f'[DONE!] {record.id}')
        
