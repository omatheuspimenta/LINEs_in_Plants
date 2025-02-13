import os
from fnmatch import fnmatch


input_folder = r'../data/LINES/'
input_domains = r'../data/domains/'


file_names = [name for name in os.listdir(input_folder) if fnmatch(name, "*.fasta")]
domains_names = [name for name in os.listdir(input_domains)]


for file in file_names:
    file_splited = file.split('.fasta')[0].split('_') 
    file_name = '_'.join(file_splited[0:-2])
    if file_name in domains_names:
        output_folder = input_domains + file_name + '/LINES/'
        os.makedirs(output_folder, exist_ok=True)
        
        input_file = input_folder + file
        output_file = output_folder + file
        # src_path = os.path.join(input_file, file)
        # dest_path = os.path.join(output_folder, file)
        if os.path.exists(input_file):
            os.rename(input_file, output_file)
