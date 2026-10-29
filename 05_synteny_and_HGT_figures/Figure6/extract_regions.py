import os
import glob
import subprocess
from Bio import SeqIO

query_fasta = "etp4.fasta"
padding = 15000

# Find all Bakta fna files in the subdirectories
fna_files = glob.glob("*/*.fna")

for fna_path in fna_files:
    dir_name = os.path.dirname(fna_path)
    genome_name = os.path.basename(dir_name)
    
    # 1. Run BLAST directly against the specific .fna file
    blast_cmd = [
        "blastn", 
        "-query", query_fasta, 
        "-subject", fna_path, 
        "-outfmt", "6 sstart send"
    ]
    
    result = subprocess.run(blast_cmd, capture_output=True, text=True)
    lines = result.stdout.strip().split("\n")
    
    if not lines or lines[0] == "":
        continue
        
    # 2. Extract minimum and maximum coordinates from all hit blocks
    coords = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) == 2:
            coords.extend([int(parts[0]), int(parts[1])])
            
    if not coords:
        continue
        
    min_pos = min(coords)
    max_pos = max(coords)
    
    # 3. Locate the corresponding Bakta .gbff file
    gbff_files = glob.glob(f"{dir_name}/*.gbff")
    if not gbff_files:
        print(f"No .gbff file found in {dir_name}")
        continue
        
    gbff_path = gbff_files[0]
    
    # Calculate slice boundaries, ensuring the start is not negative
    start_slice = max(0, min_pos - padding)
    end_slice = max_pos + padding
    
    print(f"Processing {genome_name}...")
    print(f"Prophage coordinates: {min_pos} - {max_pos}")
    print(f"Slicing region: {start_slice} - {end_slice}")
    
    # 4. Read the full GenBank file, slice, and save
    record = SeqIO.read(gbff_path, "genbank")
    
    # Ensure the end boundary does not exceed the total contig length
    end_slice = min(end_slice, len(record.seq))
    slice_record = record[start_slice:end_slice]
    
    out_name = f"{genome_name}_etp4_region.gbk"
    SeqIO.write(slice_record, out_name, "genbank")
    print(f"Saved to {out_name}\n")
