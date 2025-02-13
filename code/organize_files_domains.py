import os
from fnmatch import fnmatch
import re
from collections import defaultdict

# Input and output paths
input_folder1 = r'../data/test/'
input_folder_names = [name for name in os.listdir(input_folder1)]

for inputs in input_folder_names:
    input_folder = input_folder1 + inputs
    output_folder = input_folder + '/domains/'  
    os.makedirs(output_folder, exist_ok=True)
    
    file_names = [name for name in os.listdir(input_folder) if fnmatch(name, "*.fasta")]
    
    
    # A dictionary to group filenames by domain
    domain_files = defaultdict(list)
    
    # Regex pattern to extract the domain from the filename
    pattern = re.compile(r"^(.*?)(Monocot|Bryophyta|Eudicot|Pinophyta|Polypodiophyta|Hepatophyta)_")
    
    # Group files by domain
    for file_name in file_names:
        match = pattern.match(file_name)
        if match:
            domain = match.group(1).strip().replace(" ", "_")  # Normalize domain name
            domain_files[domain].append(file_name)
    
    temp_dict = defaultdict(list)
    for domain, files in domain_files.items():
        if 'EEP' in domain:
            temp_dict['EEP'].extend(files)
        elif 'RH' in domain:
            temp_dict['RH'].extend(files)
        elif 'RT' in domain:
            temp_dict['RT'].extend(files)
        elif 'GTT' in domain:
            temp_dict['(GTT)n'].extend(files)
        elif 'DUF' in domain:
            temp_dict['ORF1_DUF'].extend(files)
        elif 'RRM' in domain:
            temp_dict['ORF1_RRM'].extend(files)
        elif 'PHA' in domain:
            temp_dict['ORF1_PHA'].extend(files)
        elif 'Zf' in domain:
            temp_dict['ORF1_Zf'].extend(files)
        elif 'Poli' in domain:
            temp_dict['PoliA'].extend(files)
        elif 'ZF' in domain:
            temp_dict['ZF'].extend(files)
    
    
    # Combine files for each domain
    for domain, files in temp_dict.items():
        combined_file_path = os.path.join(output_folder, f"{domain}.fasta")
        with open(combined_file_path, "w") as combined_file:
            for file_name in files:
                file_path = os.path.join(input_folder, file_name)
                if os.path.exists(file_path):
                    with open(file_path, "r") as f:
                        combined_file.write(f.read())
    
    print("Files have been grouped and combined into domain-specific .fasta files.")
