# High-throughput DNA engineering by mating bacteria
Takeshi Matsui*, Po-Hsiang Hung*, Han Mei, Xianan Liu, Fangfei Li, John Collins, Weiyi Li, Darach Miller, Neil Wilson, Esteban Toro, Geoffrey J. Taghon, Gavin Sherlock, Sasha Levy#  

*Equal Contribution

#Correspondence to sasha@bacstitchdna.com

## Abstract

To reduce the operational friction and scale DNA engineering, we report here an in vivo DNA assembly technology platform called SCRIVENER (S​equential ​C​onjugation and ​R​ecombination for I​n ​V​ivo ​E​longation of ​N​ucleotides with low ​ER​rors). SCRIVENER combines bacterial conjugation, in vivo DNA cutting, and in vivo homologous recombination to seamlessly stitch blocks of DNA together by mating E. coli in large arrays or pools. This workflow is simpler, cheaper, and higher throughput than current DNA assembly approaches that require DNA to be moved in and out of cells at different procedural steps. We perform over 5,000 assemblies with two to 13 DNA blocks that range from 240 bp to 8 kb and show that SCRIVENER is capable of assembling constructs as long as 23 kb at relatively high throughput and fidelity. Most SCRIVENER errors are deletions between long interspersed repeats. However, SCRIVENER can overcome these errors by enabling assembly and sequence verification at high replication at a nominal additional cost per replicate. We show that SCRIVENER can be used to build combinatorial libraries in arrays or pools, and that DNA blocks onboarded into the platform can be repurposed and reused with any other DNA block in high throughput without a PCR step. Because of these features, DNA engineering with SCRIVENER has the potential to accelerate design-build-test-learn cycles of DNA products. 


## Table of contents
1. ### Nanoplex

	Nanoplex is a specialized Python-based pipeline for nanopore sequencing analysis. It processes long-read plasmid DNA sequencing data from Oxford Nanopore Technologies (ONT) through a series of steps: custom demultiplexing, reference mapping, variant calling, and comprehensive QC assessments. This includes generating QC reports for any mutations or alterations detected in the consensus sequence (vs. the expected reference) formed by aggregating nanopore reads. Additionally, it involves dynamic clustering analysis of the lengths of the input FASTQ reads and provides quantitative assessments of sample purity. 

	The pipeline uses subprocess calls and multiprocessing to conduct parallel per-read analyses, accommodating the demands of large, highly-multiplexed plasmid DNA datasets effectively. It features a command-line interface (CLI) with customizable options, extensive inline documentation, and built-in logging to enhance the transparency and comprehensibility of the bioinformatics processes involved.

2. ### Deletion_junction_BLAST

	Include scripts to assess whether deletions that occured during DNA assembly are more likely to occur in regions with repetitive sequences.

3. ### Pooled_GPCR_analysis

   	Include scripts to align long-read plasmid DNA sequencing data from Oxford Nanopore to a set of reference sequences. Used to analyze how many of each GPCR designs are present following pooled GPCR assembly.


## How to use Nanoplex

- [Pipeline Steps](#pipeline-steps)
- [Usage](#usage)
- [Installation](#installation)
- [Dependencies](#dependencies)

## Pipeline Steps

> Note: Input is assumed to be `dorado`-basecalled ONT native barcode-demultiplexed nanopore `*.fastq[.gz]` files.

### 1. Primary - Custom Demultiplexing

Separate reads based on custom barcodes using a tailored, efficient demultiplexing program.

### 2. Secondary - Reference Mapping & Variant Calling

Align reads by mapping with `minimap2` and call variants using:
 - `bcftools` (for small mutants &/or indels)
 - `sniffles` (for large structural variants)

Core bioinformatics pipeline of commands currently being used for the secondary analysis:
```bash
minimap2 -ax map-ont -o {fastq}.sam {ref}.fasta {fastq}.fastq
samtools sort -O bam -o {fastq}.bam -T /tmp/samtools_sorting {fastq}.sam
samtools addreplacerg -r ID:{fastq} -r SM:{fastq} -o {fastq}_RGID.bam {fastq}.bam
samtools index {fastq}.bam
bcftools mpileup -Ob -Q1 --threads 8 -d 10000 -X ont --indels-cns -B --max-BQ 35 \
                 --delta-BQ 99 -F0.2 -o15 -e1 -h110 --del-bias 0.4 --indel-bias 0.7 \
                 --poly-mqual --seqq-offset 130 --indel-size 80 \
                 -a "FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/NMBZ,FORMAT/QS,"  \ 
                    "FORMAT/SP,FORMAT/SCR,INFO/AD,INFO/ADF,INFO/ADR,INFO/BQBZ,INFO/FS," \
                    "INFO/IDV,INFO/IMF,INFO/MIN_PL_SUM,INFO/MQ0F,INFO/MQBZ,INFO/NM,"    \
                    "INFO/NMBZ,INFO/RPBZ,INFO/SCBZ,INFO/SCR,INFO/SGB,INFO/VDB"
                 -o {fastq}.bcf --fasta-ref {fastq}.fasta {fastq}.bam
bcftools call -v -m --ploidy 1 --threads 8 -Oz -o {fastq}.vcf.gz {fastq}.bcf
sniffles -i {fastq}.bam -v {fastq}.vcf
gunzip {fastq}.vcf.gz
rm {fastq}.bam # 'Read Group' IDs & Sample Names were added to {fastq}_RGID.bam
```

### 3. Tertiary - More Advanced Analyses & Final QC

- Clustering and partitioning of reads by length (followed by re-running secondary analysis)
- Purity estimation calculated from read length peak-isolated population statistics
- Microhomology-Mediated Eng Joining (MMEJ) analysis to identify large deletions (not used in this study)


## Usage

### Command-Line Interface (CLI)
```bash
:/nanoplex$ python nanoplex.py -h
usage: nanoplex.py [-h] -i INPUT_FASTQ -b BARCODE_FASTA -n RUN_NAME [-p PLASMID_REF] [-j] [-v] [-y] [-d | -s]

Ultra-fast, high-accuracy automated nanopore QC pipeline.

optional arguments:
  -h, --help            show this help message and exit
  -p PLASMID_REF, --plasmid-ref PLASMID_REF
                        path to csv file containing sufficient data for mapping references to samples
  -j, --no-junctions    skip the deletions junctions auto-detection (& MMEJ) analysis
  -v, --debug           enable debug (+'verbose') mode
  -y, --save-output     save output results; do not prompt asking
  -d, --demux-only      only run the demultiplexing analysis
  -s, --skip-demux      begin with the secondary analysis (i.e., skip demultiplexing)

required arguments:
  -i INPUT_FASTQ, --input-fastq INPUT_FASTQ
                        path to the input fastq[.gz] directory
  -b BARCODE_FASTA, --barcode-fasta BARCODE_FASTA
                        path to the barcode sequences in fasta format
  -n RUN_NAME, --run-name RUN_NAME
                        name to append to output results directory

🧬 Ⓝ ᴬᴺ🦠Pˡᵉ᥊ׁׅ ・| BacStitch DNA, Inc. © 2024
```

### Run the pipeline:

```bash
python nanoplex.py \
       -i [path to directory of input.fastq] \
       -b [path to barcodes.fasta] \
       -p [path to plasmid_references_info.csv] \
       -n [run-specific name (e.g., 'chrI' or 'BGC_assembly_13')]
```


### Example Usage
Below is an example of the output directly printed to the terminal when running Nanoplex. Example input files can be found under Nanoplex/data/example_input

> Note: Run the pipeline with the `-v` option to trigger comprehensive logging output.

```bash
nanoplex$ python nanoplex.py -y \
		-i /var/lib/minknow/data/20240328_BGC_TWIST_fluroscent_Chr_linker/20240328_BGC_TWIST_fluroscent_Chr_linker/20240328_1440_MN24328_FAY53979_683cac77/fastq_pass/barcode10 \
		-b ./data/barcodes/eaBSD3_BC.fasta \
		-p /scratch/nanoplex/data/scrivener/20240328_BGC_TWIST_fluroscent_Chr_linker_refs.csv \
		-n chrI
2024-05-29 12:30:30,770     INFO [nanoplex.py:521:run_pipeline] 
 ◖𓇃𓇃𓇃𓇃˥𔖼⎡𓇃𓇃࿅𓇃𓊢𓇃𓇃𓇃𓇃ⴶ〰⸅|ꗷꗺꗷ|⸄〰ж𓇃𓇃𓇃𓇃𓇃𓏟𔘥𔖼ꗺ𔖼𔘥𓏞𓇃𓇃𓇃𔔏𓇃𓆬𔒻𔖼𔖼𓋏𔖼𔖼𔒻𓆬𓇃𓇃𓇃𓇊⸡𓇃𓇃𓇃𓇃𓉽𓇃ண⎤ꗷ𔖿ꗷ꜒𓇃𓇃𓇃𓇃𓇃◗ 
                                                                                              
                                      ⌁🧬 Ⓝ ᴬᴺ🦠Pˡᵉ᥊ׁׅ  ・                                       
                                                                                              
                            Nanopore Sequencing Analysis Pipeline                             
                          (RUNID=NPX202405291230-v1.2.22-Alpha_chrI)                          
                                                                                              
                                     Version 1.2.22-Alpha                                     
                                  BacStitch DNA, Inc. © 2024                                  
                                                                                              
◖𓇃𓇃𓇃⸠⎫ꗺ⎧⸡𓇃𓇃𓇃𓇃⸣⸠𔖼𔖼⸡⸢𓇃𓇃𓇃𓍰𓇃𓇃𓇃⎨ꗷ𔖿ꗷ⎬𓇃𓇃𓇃𓇃𓇃𓇃𓇃𓋳𓇃𓇃𓇃𓇃𓇃╗𔖼𔖼╔𓇃𓇃𓇃𓇃𓇃ཀꗷꗺꗷཫ𓇃𓇃𓇃𓇃𓇃⦄༽⸶𔖿ꗺ𔖿⸷༼⦃𓇃𓇃𓇃𓇃◗
2024-05-29 12:30:30,770     INFO [nanoplex.py:522:run_pipeline] Nanoplex analysis initiated - Version: 1.2.22-Alpha
2024-05-29 12:30:30,770   NOTICE [nanoplex.py:523:run_pipeline] ➜ Output dir = ./results/NPX202405291230-v1.2.22-Alpha_chrI
2024-05-29 12:30:30,770   NOTICE [nanoplex.py:524:run_pipeline] ➜ Log file = ./results/NPX202405291230-v1.2.22-Alpha_chrI/Pipeline_NPX202405291230-v1.2.22-Alpha_chrI.log
2024-05-29 12:30:49,898     INFO [__init__.py:579:read_and_merge_fastq] 37 FASTQ files merged from input directory.
2024-05-29 12:30:51,769     INFO [nanoplex.py:590:run_pipeline] ➲ STEP 1: PRIMARY ANALYSIS - DEMULTIPLEXING
2024-05-29 12:31:12,847     INFO [primary.py:265:demultiplex_fastq] Demultiplexing reads (via 96 barcodes)...
2024-05-29 12:31:14,524   NOTICE [primary.py:282:demultiplex_fastq] 24331 (out of 146679) sequences matched a barcode (16.588 %)
2024-05-29 12:31:14,648  SUCCESS [nanoplex.py:610:run_pipeline] ➟ ✅ SUCCESS (Step #1/3 - PRIMARY)
2024-05-29 12:31:14,648     INFO [nanoplex.py:614:run_pipeline] ➲ STEP 2: SECONDARY ANALYSIS - MAIN PIPELINE (REF MAPPING, VARIANT CALLING)
2024-05-29 12:31:17,575  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_06_GTTGGAATCCAAACCAAACC
2024-05-29 12:31:17,576  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_E_10_ACGAATTACATTAGCTAAAG
2024-05-29 12:31:17,580  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_C_08_GCATAGCCCAACTAAGACGG
2024-05-29 12:31:17,582  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_05_AGTGAACTCCTAATCGACAC
2024-05-29 12:31:17,583  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_E_01_CGTTCAGTTTTGCAAATTAT
2024-05-29 12:31:17,584  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_01_GGTATTCCACTCAGCAAAGA
2024-05-29 12:31:17,584  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_03_AATATCCAATGTTGCTCCGA
2024-05-29 12:31:17,585  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_E_11_TCACTGTTTGCGACTATTTC
2024-05-29 12:31:17,586  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_C_11_ATCCCCAATGCCCGCACGAA
2024-05-29 12:31:17,585  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_C_04_CTTCCAGTCGACTTACAAGG
2024-05-29 12:31:17,586  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_C_10_CACTGAGAGCGTAAGATTAT
2024-05-29 12:31:17,587  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_E_09_CCAGAATGCACTAATTTAAA
2024-05-29 12:31:17,587  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_C_03_GGTCACCTCGGCAGTCAACT
2024-05-29 12:31:17,589  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_E_08_TACTAAAACAAGCTGATGAA
2024-05-29 12:31:17,589  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_C_01_TAAGTCTTAGGGATTTAAGA
2024-05-29 12:31:17,590  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_C_02_ACCACGCACCAAAAATGATA
2024-05-29 12:31:17,590  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_E_04_TACAGTAAAAGAGAGCCTTT
2024-05-29 12:31:17,593  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_E_06_AAAACAAAGTAGGACTGGAC
2024-05-29 12:31:17,594  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_E_05_TGAAGTCCTTACAGGATGCT
2024-05-29 12:31:17,599  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_02_TCCCTCGGTAGATTGATGTT
2024-05-29 12:31:17,692  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_C_05_AATTATGCAGAAACTTAAGA
2024-05-29 12:31:17,693  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_C_12_CCACCTTTTCAAGAGTTCAT
2024-05-29 12:31:17,693  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_E_03_CTTGATGAAGTGGTTACGAC
2024-05-29 12:31:17,693  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_E_12_CACATGTGTGCCCCATCAAC
2024-05-29 12:31:17,698  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_C_07_ATACCGCTAGATATAATTTA
2024-05-29 12:31:30,335  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_01_GGTATTCCACTCAGCAAAGA
2024-05-29 12:31:30,365  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_05_AGTGAACTCCTAATCGACAC
2024-05-29 12:31:30,766  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_03_AATATCCAATGTTGCTCCGA
2024-05-29 12:31:30,769  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_06_GTTGGAATCCAAACCAAACC
2024-05-29 12:31:30,878  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_02_TCCCTCGGTAGATTGATGTT
2024-05-29 12:31:40,243  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_01_ACACATGAAAATGACCGACG
2024-05-29 12:31:40,709  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_12_GCTCCTGAATCGGCGAAAAT
2024-05-29 12:31:41,612  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_03_CGAAGTCATCTCAAGATGTA
2024-05-29 12:31:42,187  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_04_GTACAGATTTGATCGAATAA
2024-05-29 12:31:42,967  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_05_AGGACGGAAAAGGGCCTTGT
2024-05-29 12:31:43,125  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_04_AGGGATCTACGGCTTAAGCT
2024-05-29 12:31:44,160  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_12_ATTCAGGGTTAGAAGACAGA
2024-05-29 12:31:44,260  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_06_CACCTCGTTGCCAATGCCAT
2024-05-29 12:31:44,383  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_11_AATTAGTAGAGGAGGAGACG
2024-05-29 12:31:44,560  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_07_TATCAATTCAAACATGGAAC
2024-05-29 12:31:44,874  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_02_AAAACCATGTATATGCGAGG
2024-05-29 12:31:45,168  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_09_GACTTTCCTCAGGACCTATT
2024-05-29 12:31:45,295  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_10_TAGACTTCATGCTGCGTATC
2024-05-29 12:31:45,302  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_09_CAAACACCAGTAGCATCCGT
2024-05-29 12:31:46,254  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_11_ACAAAGCATGTGGAAGCCCA
2024-05-29 12:31:46,580  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_08_TATCACAGACCCTGACAATA
2024-05-29 12:31:47,123  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_10_TATGGAATATCTGAATACCA
2024-05-29 12:31:49,174  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_07_ACTATGAAGCCCATATGATA
2024-05-29 12:31:50,525  SUCCESS [secondary.py:435:process_fastq_sample] ⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_08_GTTTATCGCGTAACGATACT
2024-05-29 12:31:50,841   NOTICE [nanoplex.py:627:run_pipeline] Note: 72 empty FASTQ files.
2024-05-29 12:31:50,841  SUCCESS [nanoplex.py:630:run_pipeline] ➟ ✅ SUCCESS (Step #2/3 - SECONDARY)
2024-05-29 12:31:50,842     INFO [nanoplex.py:634:run_pipeline] ➲ STEP 3: TERTIARY ANALYSIS - CLUSTERING, PURITY, QC SUMMARY
2024-05-29 12:31:50,842     INFO [nanoplex.py:640:run_pipeline] Performing read lengths clustering & purity estimation analysis...
2024-05-29 12:32:01,279  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_02_TCCCTCGGTAGATTGATGTT
2024-05-29 12:32:01,563  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_05_AGTGAACTCCTAATCGACAC
2024-05-29 12:32:01,616  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_03_AATATCCAATGTTGCTCCGA
2024-05-29 12:32:02,167  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_03_CGAAGTCATCTCAAGATGTA
2024-05-29 12:32:02,463  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_06_CACCTCGTTGCCAATGCCAT
2024-05-29 12:32:02,616  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_12_ATTCAGGGTTAGAAGACAGA
2024-05-29 12:32:03,382  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_09_CAAACACCAGTAGCATCCGT
2024-05-29 12:32:03,871  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_01_GGTATTCCACTCAGCAAAGA
2024-05-29 12:32:04,081  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_06_GTTGGAATCCAAACCAAACC
2024-05-29 12:32:04,611  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_02_AAAACCATGTATATGCGAGG
2024-05-29 12:32:04,985  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_01_ACACATGAAAATGACCGACG
2024-05-29 12:32:05,242  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_05_AGGACGGAAAAGGGCCTTGT
2024-05-29 12:32:05,305  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_12_GCTCCTGAATCGGCGAAAAT
2024-05-29 12:32:05,337  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_04_AGGGATCTACGGCTTAAGCT
2024-05-29 12:32:05,390  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_09_GACTTTCCTCAGGACCTATT
2024-05-29 12:32:05,999  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_04_GTACAGATTTGATCGAATAA
2024-05-29 12:32:06,262  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_07_TATCAATTCAAACATGGAAC
2024-05-29 12:32:06,741  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_11_ACAAAGCATGTGGAAGCCCA
2024-05-29 12:32:06,859  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_11_AATTAGTAGAGGAGGAGACG
2024-05-29 12:32:07,535  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_10_TAGACTTCATGCTGCGTATC
2024-05-29 12:32:07,583  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_08_TATCACAGACCCTGACAATA
2024-05-29 12:32:08,354  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_10_TATGGAATATCTGAATACCA
2024-05-29 12:32:13,301  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_07_ACTATGAAGCCCATATGATA
2024-05-29 12:32:13,964  SUCCESS [tertiary.py:961:analyze_mmej_and_plot_coverage] ⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_08_GTTTATCGCGTAACGATACT
2024-05-29 12:32:14,334     INFO [nanoplex.py:690:run_pipeline] RE-RUNNING SECONDARY ANALYSIS on clustered fastq files...
2024-05-29 12:32:17,250  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_A_12_GCTCCTGAATCGGCGAAAAT_cluster_1
2024-05-29 12:32:17,250  WARNING [secondary.py:135:process_fastq_sample] ⚠︎ <100 total reads for: FAY53979_pass_barcode10_683cac77_3c60b55a_merged_eaBSD3_B_05_AGGACGGAAAAGGGCCTTGT_cluster_1
2024-05-29 12:32:28,414   NOTICE [nanoplex.py:721:run_pipeline] Mutations/indels detected in samples: 

                          A05                                                                            

2024-05-29 12:32:28,465  SUCCESS [nanoplex.py:834:run_pipeline] ➟ ✅ SUCCESS (Step #3/3 - TERTIARY)
2024-05-29 12:32:28,465  SUCCESS [nanoplex.py:844:run_pipeline] 🏁 Nanoplex v1.2.22-Alpha «NPX202405291230-v1.2.22-Alpha_chrI» analysis complete. 

	⌛️ Total execution time: 1 min, 57.694 s 

```


## Installation

Clone the repository and ensure dependencies are installed (see below).

1. Clone the repository

2. Install dependencies (see [Dependencies](#dependencies)).

```bash
cd nanoplex/
./install_tools.sh
```

## Dependencies

Nanoplex requires the following dependencies:

- Python 3.6+
- samtools
- bcftools
- minimap2
- sniffles

Install these dependencies using appropriate package managers. (See the `./install_tools.sh` script.)

## Raw sequence data

Raw sequence data for this study can be found on Sequence Read Archive (SRA BioProject PRJNA1150152) 
