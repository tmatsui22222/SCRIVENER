####read the supplemental data with the expected plasmid sequence of each construct
data = read.table('~/Desktop/SCRIVENER_paper/Supplemental_data/Supplementary_Data_3.tsv', header = T, sep = '\t')

###read the supplmentary table with all the variants identified
muts = read.table('~/Desktop/SCRIVENER_paper/Supplemental_data/Table_S5.tsv', header = T, sep = '\t')

###filter data to only include full length assemblies. Does not include any deletions identified in BGC partial assemblies, GPCR designs, or bridge assemblies
bgc5 = data[grepl('stitch_number5', data$Construct),]
bgc9 = data[grepl('stitch_number9', data$Construct),]
bgc11 = data[grepl('stitch_number11', data$Construct),]
bgc13 = data[grepl('stitch_number13', data$Construct),]
yeast = data[grepl('Yeast', data$Construct),]
syn = data[grepl('Synthetic', data$Construct),]
full = rbind(bgc13, yeast, syn, bgc9[grepl('Neocamarosporium', bgc9$Construct),], bgc5[grepl('Cyanothece', bgc5$Construct),], bgc5[grepl('Scytonema', bgc5$Construct),], bgc11[grepl('flavus', bgc11$Construct),])

###isolate only deletions identified the full length assemblies.
muts = muts[which(muts$Construct %in% full$Construct),]
del = muts[which(muts$variant_type == 'deletion'),]
del$info = paste(del$well, del$start, del$end, del$reference, del$alt, sep = '_')
del = del[!duplicated(del$info),]

###extract the 101bp windows surrounding the start and end positions of each deletion identified. For each deletion, randomly sample 100 random pairs of 101 bp windows, separated by the same distance as the deletion, from the construct where the deletion was detected.
targets2 = sort(unique(del$Construct))

data = data.frame()
for(i in 1:nrow(del)){
	temp = del[i,]
	file = full[which(full$Construct == temp$Construct), 'Expected_plasmid_sequence']
	start = (temp$start - 50):(temp$start + 50)
	end = (temp$end - 50):(temp$end + 50)
	blastl = substr(file, min(start), max(start))
	blastr = substr(file, min(end), max(end))
	data = rbind(data, c(temp$Construct, blastl, blastr))
	
	char = nchar(file)
	for(q in 1:100){
		s = sample(char-abs(temp$variant_length)-51, 1)
		e = s + abs(temp$variant_length)
		start = (s - 50):(s + 50)
		end = (e - 50):(e + 50)
		
		blastl = substr(file, min(start), max(start))
		blastr = substr(file, min(end), max(end))
		data = rbind(data, c(paste0(temp$Construct, '_p', q), blastl, blastr))
	}
	
}
names(data) = c('Target', 'blast_left', 'blast_right')
write.table(data, file = '~/Desktop/SCRIVENER_paper/Supplemental_data/only_full_assembly_blast.tsv', sep = '\t', row.names = F, quote = F)

#########################################################################################################################################################
##run deletion_junction_blast.sh to get the max BLAST score for each pair of 101bp windows 
##requires BLAST to be installed (https://www.ncbi.nlm.nih.gov/books/NBK569861/)
#########################################################################################################################################################

###read the output table after running deletion_junction_blast.sh
table = read.table('~/Desktop/SCRIVENER_paper/Supplemental_data/only_full_assembly_blast_out.tsv', header = T, as.is = T, sep = '\t')
###separate data based on actual deletion junctions from permuted junctions
sample = table[seq(1, nrow(table), 101),]
perm = table[-seq(1, nrow(table), 101),]

sample2 = perm$Max_BLAST_Score
sample1 = sample$Max_BLAST_Score
###one-sided wilcoxon rank sum test to test if actual deletion junctions have a larger BLAST scores than random permuted junctions
pval = wilcox.test(sample1,sample2, alternative = 'greater')$p.value





