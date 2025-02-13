#!/bin/bash

# Find all .fasta files and compress each into its own .tar.xz archive
for file in *.fasta; do
    if [ -f "$file" ]; then
        tar -cJf "${file}.tar.xz" "$file"
        echo "Compressed: ${file} -> ${file}.tar.xz"
    fi
done

echo "Compression process completed."

rm *.fasta

echo "Deleted all fasta files."
