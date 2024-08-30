#!/bin/bash

stitch_orderL=("Stitch1_replicate" "Stitch2_replicate" "Stitch3_replicate" "Stitch4_replicate")
repeatL=("1-1" "1-2" "2-1" "2-2")
refL=("All_6_GPCR_combinations_1st_stitch.fasta" "All_6_GPCR_combinations_2nd_stitch.fasta" "All_6_GPCR_combinations_3rd_stitch.fasta" "All_6_GPCR_combinations_4th_stitch.fasta")

# 1st stitch
for i in {0..3}; do
    for repeat in "${repeatL[@]}"; do
        mkdir -p ./raw_demuxed_file/merged_reads/sam
        mkdir -p ./raw_demuxed_file/merged_reads/paf
        minimap2 -ax map-ont "${refL[i]}" "./raw_demuxed_file/merged_reads/${stitch_orderL[i]}${repeat}.fastq" > "./raw_demuxed_file/merged_reads/sam/${stitch_orderL[i]}${repeat}.sam"
        minimap2 -c "${refL[i]}" "./raw_demuxed_file/merged_reads/${stitch_orderL[i]}${repeat}.fastq" > "./raw_demuxed_file/merged_reads/paf/${stitch_orderL[i]}${repeat}.paf"
    done
done
        
#run in commed by
#chmod +x pooled_GPCR_run_minimap2.sh
#./pooled_GPCR_run_minimap2.sh