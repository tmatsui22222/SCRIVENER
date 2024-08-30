#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.tsv output.tsv"
    exit 1
fi

input_file=$1
output_file=$2

# Create output file and write the header
echo -e "Name\tMax_BLAST_Score\tBest_Match_Sequence\tLength\tPercent_Identity" > "$output_file"

# Get the total number of lines in the input file (excluding the header)
total_lines=$(tail -n +2 "$input_file" | wc -l)
current_line=0

# Function to update the progress bar
update_progress() {
    current_line=$((current_line + 1))
    progress=$(echo "scale=2; $current_line / $total_lines * 100" | bc)
    printf "\rProgress: %d/%d (%.2f%%)" "$current_line" "$total_lines" "$progress"
}

# Function to process each line of the input file
process_line() {
    name=$1
    seq1=$2
    seq2=$3

    # Create unique temporary files for the sequences in /var/tmp
    seq1_file=$(mktemp /var/tmp/seq1_XXXXXX.fasta) || { echo "Failed to create temporary file for $name"; return 1; }
    seq2_file=$(mktemp /var/tmp/seq2_XXXXXX.fasta) || { echo "Failed to create temporary file for $name"; rm -f "$seq1_file"; return 1; }

    echo -e ">seq1\n$seq1" > "$seq1_file"
    echo -e ">seq2\n$seq2" > "$seq2_file"

    # Create a BLAST database for seq1
    makeblastdb -in "$seq1_file" -dbtype nucl -out seq1_db -logfile /dev/null
    if [ $? -ne 0 ]; then
        echo "Error creating BLAST database for $name"
        cat "$seq1_file"
        rm -f "$seq1_file" "$seq2_file" seq1_db.*
        return 1
    fi

    # Run BLAST comparing seq2 to seq1 database with parameters for short sequences
    blastn -query "$seq2_file" -db seq1_db -outfmt "6 qseqid sseqid bitscore length pident qseq sseq" -word_size 4 -reward 1 -penalty -1 -gapopen 5 -gapextend 2 -out blast_results.txt
    if [ $? -ne 0 ]; then
        echo "Error running BLAST for $name"
        cat "$seq2_file"
        rm -f "$seq1_file" "$seq2_file" seq1_db.* blast_results.txt
        return 1
    fi

    # Get the maximum BLAST score and the corresponding best match sequence
    max_score=0
    best_match_sequence=""
    best_length=0
    best_pident=0
    while IFS=$'\t' read -r qseqid sseqid bitscore length pident qseq sseq
    do
        if (( $(echo "$bitscore > $max_score" | bc -l) )); then
            max_score=$bitscore
            best_match_sequence=$sseq
            best_length=$length
            best_pident=$pident
        fi
    done < blast_results.txt

    # Write the result to the output file
    echo -e "$name\t$max_score\t$best_match_sequence\t$best_length\t$best_pident" >> "$output_file"

    # Clean up temporary files
    rm -f "$seq1_file" "$seq2_file" seq1_db.* blast_results.txt
    
    # Update progress bar
    update_progress
}

# Read the TSV file line by line and process each line
tail -n +2 "$input_file" | while IFS=$'\t' read -r name seq1 seq2; do
    process_line "$name" "$seq1" "$seq2"
done

# Print a newline after the progress bar
echo
