import os
import re
from collections import defaultdict

# Input and output folder paths
input_folder = r'../data/domains/'
output_folder = r'../data/test/'

# Ensure the output folder exists
os.makedirs(output_folder, exist_ok=True)

# Regex to extract organism and domain information
pattern = re.compile(r"(?:\w+)?(Bryophyta_[\w]+|Monocot_[\w]+|Eudicot_[\w]+|Pinophyta_[\w]+|Polypodiophyta_[\w]+|Hepatophyta_[\w]+)_(L\d+|RTE|I)_\d+\.fasta")

# Process filenames
with open("file_names.txt", "r") as file:
    file_names = file.read().splitlines()

# Organize files by organism
organism_files = defaultdict(list)

for file_name in file_names:
    match = pattern.search(file_name)
    if match:
        organism = match.group(1)  # Extract the organism name
        organism_files[organism].append(file_name)

# Move and organize files
for organism, files in organism_files.items():
    organism_dir = os.path.join(output_folder, organism)
    os.makedirs(organism_dir, exist_ok=True)

    # Move files into the organism folder
    for file_name in files:
        src_path = os.path.join(input_folder, file_name)
        dest_path = os.path.join(organism_dir, file_name)
        if os.path.exists(src_path):
            os.rename(src_path, dest_path)

print("All files have been moved into organism-specific folders.")
