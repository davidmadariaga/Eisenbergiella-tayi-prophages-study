#!/bin/bash

# =========================================================
# DIAMOND CHECK: DefenseFinder vs Prophages (Sequence Match)
# =========================================================

# 1. SETUP - EXACT PATHS
# ----------------------
# This is the path you found where the Bakta .faa files live
ORIGINAL_GENOME_PATH="/home/d/Documents/research/eisenbergiella_tayi/bakta_genomes/selectedannotatedgenomes"

# Thresholds advised by your senior
ID_THRESH=30    # 30% Identity
COV_THRESH=70   # 70% Coverage

mkdir -p diamond_analysis
OUTPUT_FILE="diamond_analysis/final_overlap_results.tsv"

# Write Output Header
echo -e "Genome\tProphage_Protein_ID\tDefense_Hit_ID\tIdentity\tCoverage\tSystem_Type" > "$OUTPUT_FILE"

echo "Starting Analysis..."

# 2. LOOP THROUGH GENOMES
# -----------------------
# We loop through the folders in 'defensesystems' to find the inputs
for genome_dir in defensesystems/GCF_*; do
    genome_name=$(basename "$genome_dir")
    
    # A. Find the Prophage File (Query) inside the folder
    prophage_fasta=$(find "$genome_dir" -name "*provirus*.faa" | head -n 1)

    # B. Find the Original Bakta Genome (Subject)
    # We look for the .faa file in the path you found with grep
    full_genome_fasta="${ORIGINAL_GENOME_PATH}/${genome_name}/${genome_name}.faa"

    # Check if both files exist
    if [[ -f "$prophage_fasta" && -f "$full_genome_fasta" ]]; then
        echo "Processing $genome_name..."

        # C. Get the list of Defense IDs for this genome
        # We extract column 2 (Hit_id) from the DefenseFinder genes file
        df_genes_file="$genome_dir/${genome_name}_defense_finder_genes.tsv"
        
        if [[ -f "$df_genes_file" ]]; then
            tail -n +2 "$df_genes_file" | cut -f 2 > "temp_defense_ids.txt"

            # D. Extract the sequences of ONLY the defense genes
            # This creates a mini-database of defense proteins
            seqkit grep -f "temp_defense_ids.txt" "$full_genome_fasta" > "temp_defense_seqs.faa" 2>/dev/null

            # E. Run DIAMOND (BLASTP)
            # Query = Prophage Proteins
            # Subject = Defense Proteins
            diamond makedb --in "temp_defense_seqs.faa" -d "temp_db" --quiet
            
            diamond blastp \
                -q "$prophage_fasta" \
                -d "temp_db" \
                -o "temp_matches.tsv" \
                --outfmt 6 qseqid sseqid pident qcovhsp \
                --quiet

            # F. Filter Results (Identity >= 30, Coverage >= 70)
            awk -v genome="$genome_name" -v id="$ID_THRESH" -v cov="$COV_THRESH" \
            '$3 >= id && $4 >= cov {print genome "\t" $1 "\t" $2 "\t" $3 "\t" $4}' \
            "temp_matches.tsv" >> "$OUTPUT_FILE"

            rm temp_defense_ids.txt temp_defense_seqs.faa temp_db.dmnd temp_matches.tsv
        else
            echo "  Warning: DefenseFinder genes file missing for $genome_name"
        fi
    else
        echo "  Skipping $genome_name (Missing Prophage or Bakta FASTA file)"
    fi
done

echo "=== DONE ==="
echo "Results saved to: $OUTPUT_FILE"
cat "$OUTPUT_FILE"
