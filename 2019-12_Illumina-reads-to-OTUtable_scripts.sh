###This script is used to process rDNA amplicon data and generate OTU table
###Author: Chitong Rao, PhD, 2019/12, Boston Childrens's Hospital / Harvard Medical School

##############################################################################

###Demultiplex the Illumina i5/i7-demultiplexed subpools into single samples using skewer with the first 12nt exact match
##prerequisite: download and install skewer (https://github.com/relipmoc/skewer)

#in python2 environment
source activate python2
mkdir temp/ demultiplexed/ demultiplex_log/

for i in $(cut -f 1 subpool_bac16S_list | sort | uniq); do
./skewer -x ./bac16S_forward-barcode.fa -y ./bac16S_reverse-barcode.fa -M ./matrix_8F_12R --mode ap -r 0 -t 2 -b --barcode raw/$i\_*_R1_001.fastq.gz raw/$i\_*_R2_001.fastq.gz -o temp/$i
grep $i subpool_bac16S_list > subpool_temp
while IFS=$'\t' read -r -a array
do
sample=${array[1]}
code=${array[2]}
mv temp/$i\-assigned-$code\-pair1.fastq demultiplexed/$sample\_bac16S_R1.fastq
mv temp/$i\-assigned-$code\-pair2.fastq demultiplexed/$sample\_bac16S_R2.fastq
done < subpool_temp
mv temp/$i\-assigned.log demultiplex_log/$i\_demultiplex.log
rm subpool_temp temp/*
done

for i in $(cut -f 1 subpool_arch16S_list | sort | uniq); do
./skewer -x ./arch16S_forward-barcode.fa -y ./arch16S_reverse-barcode.fa -M ./matrix_8F_12R --mode ap -r 0 -t 2 -b --barcode raw/$i\_*_R1_001.fastq.gz raw/$i\_*_R2_001.fastq.gz -o temp/$i
grep $i subpool_arch16S_list > subpool_temp
while IFS=$'\t' read -r -a array
do
sample=${array[1]}
code=${array[2]}
mv temp/$i\-assigned-$code\-pair1.fastq demultiplexed/$sample\_arch16S_R1.fastq
mv temp/$i\-assigned-$code\-pair2.fastq demultiplexed/$sample\_arch16S_R2.fastq
done < subpool_temp
mv temp/$i\-assigned.log demultiplex_log/$i\_demultiplex.log
rm subpool_temp temp/*
done

for i in $(cut -f 1 subpool_ITS1_list | sort | uniq); do
./skewer -x ./ITS1_forward-barcode.fa -y ./ITS1_reverse-barcode.fa -M ./matrix_8F_12R --mode ap -r 0 -t 2 -b --barcode raw/$i\_*_R1_001.fastq.gz raw/$i\_*_R2_001.fastq.gz -o temp/$i
grep $i subpool_ITS1_list > subpool_temp
while IFS=$'\t' read -r -a array
do
sample=${array[1]}
code=${array[2]}
mv temp/$i\-assigned-$code\-pair1.fastq demultiplexed/$sample\_ITS1_R1.fastq
mv temp/$i\-assigned-$code\-pair2.fastq demultiplexed/$sample\_ITS1_R2.fastq
done < subpool_temp
mv temp/$i\-assigned.log demultiplex_log/$i\_demultiplex.log
rm subpool_temp temp/*
done

##############################################################################

###OTU table pipeline starting from demultiplexed reads
##prerequisite: download the FAST package (https://github.com/ZeweiSong/FAST) and install stand-alone cutadapt, vsearch and usearch packages in the FAST/ folder
##download appropriate rDNA database and use qiime2 (and vsearch) for OTU assignment

###Use Nextseq 200nt read 1 only (higher quality) and trim uniformly to 100nt

#bac16SV4
cp demultiplexed/*_bac16S_R1.fastq bac16S_R1/
FAST/fast.py -generate_mapping -i bac16S_R1/ -o bac16S_R1/read1_map.txt
FAST/fast.py -add_labels -m bac16S_R1/read1_map.txt -i bac16S_R1/ -o bac16S_R1/read1_labeled -t 4
FAST/fast.py -merge_seqs -i bac16S_R1/read1_labeled -o bac16S_R1/NICU_bac16S_read1.fastq
rm -r bac16S_R1/read1_labeled
#modify primer sequence as needed
cutadapt -g GTGYCAGCMGCCGCGGTAA -o bac16S_R1/NICU_bac16S_cut1.fastq bac16S_R1/NICU_bac16S_read1.fastq --discard-untrimmed
rm bac16S_R1/NICU_bac16S_read1.fastq
FAST/vsearch --fastq_filter bac16S_R1/NICU_bac16S_cut1.fastq --fastq_truncee 1 --fastaout bac16S_R1/NICU_bac16S_trimmed_filtered.fasta --fasta_width 0
rm bac16S_R1/NICU_bac16S_cut1.fastq
FAST/fast.py -stat_seqs -i bac16S_R1/NICU_bac16S_trimmed_filtered.fasta -o bac16S_R1/NICU_bac16S_trimmed_filtered_report.txt
FAST/fast.py -truncate_seqs -i bac16S_R1/NICU_bac16S_trimmed_filtered.fasta -fixed_length 100 -o  bac16S_R1/NICU_bac16S_trimmed_filtered_L100.fasta
FAST/fast.py -dereplicate -sizeout -i bac16S_R1/NICU_bac16S_trimmed_filtered_L100.fasta -o bac16S_R1/NICU_bac16S_trimmed_filtered_L100_derep -t 4
rm bac16S_R1/NICU_bac16S_trimmed_filtered_L100.fasta
FAST/usearch -unoise3 bac16S_R1/NICU_bac16S_trimmed_filtered_L100_derep.fasta -tabbedout bac16S_R1/NICU_bac16S_R1only_unoise-denoise.txt -unoise_alpha 5
grep "\tchfilter\tzotu" bac16S_R1/NICU_bac16S_R1only_unoise-denoise.txt | cut -f 1 | cut -d ";" -f 1 > bac16S_R1/temp-list.txt
FAST/fast.py -filter_otu_map -i bac16S_R1/NICU_bac16S_trimmed_filtered_L100_derep.txt -name_list bac16S_R1/temp-list.txt -o bac16S_R1/NICU_bac16S_R1only_zotu.txt
FAST/fast.py -pick_seqs -sizeout -i bac16S_R1/NICU_bac16S_trimmed_filtered_L100_derep.fasta -map bac16S_R1/NICU_bac16S_R1only_zotu.txt -o bac16S_R1/NICU_bac16S_R1only_zotu.fasta
FAST/fast.py -make_otu_table -qiime_map bac16S_R1/NICU_bac16S_R1only_zotu.txt -o bac16S_R1/NICU_bac16S_R1only_otu-table.txt
FAST/vsearch --usearch_global bac16S_R1/NICU_bac16S_R1only_zotu.fasta -db ../rDNA_database/silva_nr_v132_train_set.fa --userout bac16S_R1/NICU_bac16S_R1only_vsearch.txt --userfields query+target+ql+pairs+id --id 0.6
FAST/fast.py -filter_taxonomy -i bac16S_R1/NICU_bac16S_R1only_vsearch.txt -op bac16S_R1/temp.txt -match_length 0.9 -pident 97
perl -pe 's/;size=\d*\t/\t/' bac16S_R1/temp.txt > bac16S_R1/NICU_bac16S_R1only_vsearch_assigned.txt
FAST/fast.py -assign_taxonomy -keep_all -otu bac16S_R1/NICU_bac16S_R1only_otu-table.txt -tax bac16S_R1/NICU_bac16S_R1only_vsearch_assigned.txt -o bac16S_R1/NICU_bac16S_R1only_otu-table_annotated.txt -scores
rm bac16S_R1/*temp* bac16S_R1/*_bac16S_R1.fastq
perl -pe 's/;size=\d*//' bac16S_R1/NICU_bac16S_R1only_zotu.fasta > bac16S_R1/temp
mv bac16S_R1/temp bac16S_R1/NICU_bac16S_R1only_zotu.fasta

#arch16S
cp demultiplexed/*_arch16S_R1.fastq arch16S_R1/
FAST/fast.py -generate_mapping -i arch16S_R1/ -o arch16S_R1/read1_map.txt
FAST/fast.py -add_labels -m arch16S_R1/read1_map.txt -i arch16S_R1/ -o arch16S_R1/read1_labeled -t 4
FAST/fast.py -merge_seqs -i arch16S_R1/read1_labeled -o arch16S_R1/NICU_arch16S_read1.fastq
rm -r arch16S_R1/read1_labeled
#modify primer sequence as needed
cutadapt -g TGYCAGCCGCCGCGGTAAHACCVGC -e 0.2 -o arch16S_R1/NICU_arch16S_cut1.fastq arch16S_R1/NICU_arch16S_read1.fastq --discard-untrimmed
rm arch16S_R1/NICU_arch16S_read1.fastq
FAST/vsearch --fastq_filter arch16S_R1/NICU_arch16S_cut1.fastq --fastq_truncee 1 --fastaout arch16S_R1/NICU_arch16S_trimmed_filtered.fasta --fasta_width 0
rm arch16S_R1/NICU_arch16S_cut1.fastq
FAST/fast.py -stat_seqs -i arch16S_R1/NICU_arch16S_trimmed_filtered.fasta -o arch16S_R1/NICU_arch16S_trimmed_filtered_report.txt
FAST/fast.py -truncate_seqs -i arch16S_R1/NICU_arch16S_trimmed_filtered.fasta -fixed_length 100 -o  arch16S_R1/NICU_arch16S_trimmed_filtered_L100.fasta
FAST/fast.py -dereplicate -sizeout -i arch16S_R1/NICU_arch16S_trimmed_filtered_L100.fasta -o arch16S_R1/NICU_arch16S_trimmed_filtered_L100_derep -t 4
rm arch16S_R1/NICU_arch16S_trimmed_filtered_L100.fasta
FAST/usearch -unoise3 arch16S_R1/NICU_arch16S_trimmed_filtered_L100_derep.fasta -tabbedout arch16S_R1/NICU_arch16S_R1only_unoise-denoise.txt -unoise_alpha 5
grep "\tchfilter\tzotu" arch16S_R1/NICU_arch16S_R1only_unoise-denoise.txt | cut -f 1 | cut -d ";" -f 1 > arch16S_R1/temp-list.txt
FAST/fast.py -filter_otu_map -i arch16S_R1/NICU_arch16S_trimmed_filtered_L100_derep.txt -name_list arch16S_R1/temp-list.txt -o arch16S_R1/NICU_arch16S_R1only_zotu.txt
FAST/fast.py -pick_seqs -sizeout -i arch16S_R1/NICU_arch16S_trimmed_filtered_L100_derep.fasta -map arch16S_R1/NICU_arch16S_R1only_zotu.txt -o arch16S_R1/NICU_arch16S_R1only_zotu.fasta
FAST/fast.py -make_otu_table -qiime_map arch16S_R1/NICU_arch16S_R1only_zotu.txt -o arch16S_R1/NICU_arch16S_R1only_otu-table.txt
FAST/vsearch --usearch_global arch16S_R1/NICU_arch16S_R1only_zotu.fasta -db ../rDNA_database/silva_nr_v132_train_set.fa --userout arch16S_R1/NICU_arch16S_R1only_vsearch.txt --userfields query+target+ql+pairs+id --id 0.6
FAST/fast.py -filter_taxonomy -i arch16S_R1/NICU_arch16S_R1only_vsearch.txt -op arch16S_R1/temp.txt -match_length 0.9 -pident 97
perl -pe 's/;size=\d*\t/\t/' arch16S_R1/temp.txt > arch16S_R1/NICU_arch16S_R1only_vsearch_assigned.txt
FAST/fast.py -assign_taxonomy -keep_all -otu arch16S_R1/NICU_arch16S_R1only_otu-table.txt -tax arch16S_R1/NICU_arch16S_R1only_vsearch_assigned.txt -o arch16S_R1/NICU_arch16S_R1only_otu-table_annotated.txt -scores
rm arch16S_R1/*temp* arch16S_R1/*_arch16S_R1.fastq
perl -pe 's/;size=\d*//' arch16S_R1/NICU_arch16S_R1only_zotu.fasta > arch16S_R1/temp
mv arch16S_R1/temp arch16S_R1/NICU_arch16S_R1only_zotu.fasta

#ITS1
cp demultiplexed/*_ITS1_R1.fastq ITS1_R1/
FAST/fast.py -generate_mapping -i ITS1_R1/ -o ITS1_R1/read1_map.txt
FAST/fast.py -add_labels -m ITS1_R1/read1_map.txt -i ITS1_R1/ -o ITS1_R1/read1_labeled -t 4
FAST/fast.py -merge_seqs -i ITS1_R1/read1_labeled -o ITS1_R1/NICU_ITS1_read1.fastq
rm -r ITS1_R1/read1_labeled
cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -o ITS1_R1/NICU_ITS1_cut1.fastq ITS1_R1/NICU_ITS1_read1.fastq
rm ITS1_R1/NICU_ITS1_read1.fastq
cutadapt -g CTTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATT -e 0.2 -o ITS1_R1/NICU_ITS1_cut2.fastq ITS1_R1/NICU_ITS1_cut1.fastq --discard-untrimmed
rm ITS1_R1/NICU_ITS1_cut1.fastq
cutadapt -a AACTTTCAACAACGGATCTCTTGGYTCTSGCATCGATGAAGAACGCAGC -o ITS1_R1/NICU_ITS1_cut3.fastq ITS1_R1/NICU_ITS1_cut2.fastq
rm ITS1_R1/NICU_ITS1_cut2.fastq
FAST/vsearch --fastq_filter ITS1_R1/NICU_ITS1_cut3.fastq --fastq_truncee 1 --fastaout ITS1_R1/NICU_ITS1_trimmed_filtered.fasta --fasta_width 0
rm ITS1_R1/NICU_ITS1_cut3.fastq
FAST/fast.py -stat_seqs -i ITS1_R1/NICU_ITS1_trimmed_filtered.fasta -o ITS1_R1/NICU_ITS1_trimmed_filtered_report.txt
FAST/fast.py -truncate_seqs -i ITS1_R1/NICU_ITS1_trimmed_filtered.fasta -fixed_length 100 -o  ITS1_R1/NICU_ITS1_trimmed_filtered_L100.fasta
FAST/fast.py -dereplicate -sizeout -i ITS1_R1/NICU_ITS1_trimmed_filtered_L100.fasta -o ITS1_R1/NICU_ITS1_trimmed_filtered_L100_derep -t 4
rm ITS1_R1/NICU_ITS1_trimmed_filtered_L100.fasta
FAST/usearch -unoise3 ITS1_R1/NICU_ITS1_trimmed_filtered_L100_derep.fasta -tabbedout ITS1_R1/NICU_ITS1_R1only_unoise-denoise.txt -unoise_alpha 5
grep "\tchfilter\tzotu" ITS1_R1/NICU_ITS1_R1only_unoise-denoise.txt | cut -f 1 | cut -d ";" -f 1 > ITS1_R1/temp-list.txt
FAST/fast.py -filter_otu_map -i ITS1_R1/NICU_ITS1_trimmed_filtered_L100_derep.txt -name_list ITS1_R1/temp-list.txt -o ITS1_R1/NICU_ITS1_R1only_zotu.txt
FAST/fast.py -pick_seqs -sizeout -i ITS1_R1/NICU_ITS1_trimmed_filtered_L100_derep.fasta -map ITS1_R1/NICU_ITS1_R1only_zotu.txt -o ITS1_R1/NICU_ITS1_R1only_zotu.fasta
FAST/fast.py -make_otu_table -qiime_map ITS1_R1/NICU_ITS1_R1only_zotu.txt -o ITS1_R1/NICU_ITS1_R1only_otu-table.txt
FAST/vsearch --usearch_global ITS1_R1/NICU_ITS1_R1only_zotu.fasta -db ../rDNA_database/sh_general_release_dynamic_01.12.2017.fasta --userout ITS1_R1/NICU_ITS1_R1only_vsearch.txt --userfields query+target+ql+pairs+id --id 0.6
FAST/fast.py -filter_taxonomy -i ITS1_R1/NICU_ITS1_R1only_vsearch.txt -op ITS1_R1/temp.txt -match_length 0.7 -pident 90
perl -pe 's/;size=\d*\t/\t/' ITS1_R1/temp.txt > ITS1_R1/NICU_ITS1_R1only_vsearch_assigned.txt
FAST/fast.py -assign_taxonomy -keep_all -otu ITS1_R1/NICU_ITS1_R1only_otu-table.txt -tax ITS1_R1/NICU_ITS1_R1only_vsearch_assigned.txt -o ITS1_R1/NICU_ITS1_R1only_otu-table_annotated.txt -scores
rm ITS1_R1/*temp* ITS1_R1/*_ITS1_R1.fastq
perl -pe 's/;size=\d*//' ITS1_R1/NICU_ITS1_R1only_zotu.fasta > ITS1_R1/temp
mv ITS1_R1/temp ITS1_R1/NICU_ITS1_R1only_zotu.fasta

#OTU assignment using qiime2
source activate qiime2-2018.2
qiime tools import \
  --input-path bac16S_R1/NICU_bac16S_R1only_zotu.fasta \
  --output-path bac16S_R1/NICU_bac16S_R1only_sequences.qza \
  --type 'FeatureData[Sequence]'
qiime feature-classifier classify-sklearn \
  --i-classifier ../rDNA_database/silva-119-99-515-806-nb-classifier.qza \
  --i-reads bac16S_R1/NICU_bac16S_R1only_sequences.qza \
  --o-classification bac16S_R1/NICU_bac16S_R1only_sklearn-silva-taxonomy.qza
qiime metadata tabulate \
  --m-input-file bac16S_R1/NICU_bac16S_R1only_sklearn-silva-taxonomy.qza \
  --o-visualization bac16S_R1/NICU_bac16S_R1only_sklearn-silva-taxonomy.qzv

qiime tools import \
  --input-path arch16S_R1/NICU_arch16S_R1only_zotu.fasta \
  --output-path arch16S_R1/NICU_arch16S_R1only_sequences.qza \
  --type 'FeatureData[Sequence]'
qiime feature-classifier classify-sklearn \
  --i-classifier ../rDNA_database/silva-119-99-nb-classifier.qza \
  --i-reads arch16S_R1/NICU_arch16S_R1only_sequences.qza \
  --o-classification arch16S_R1/NICU_arch16S_R1only_sklearn-silva-taxonomy.qza
qiime metadata tabulate \
  --m-input-file arch16S_R1/NICU_arch16S_R1only_sklearn-silva-taxonomy.qza \
  --o-visualization arch16S_R1/NICU_arch16S_R1only_sklearn-silva-taxonomy.qzv

qiime tools import \
  --input-path ITS1_R1/NICU_ITS1_R1only_zotu.fasta \
  --output-path ITS1_R1/NICU_ITS1_R1only_sequences.qza \
  --type 'FeatureData[Sequence]'
qiime feature-classifier classify-sklearn \
  --i-classifier ../rDNA_database/UNITE_s_01.12.2017_classifier.qza \
  --i-reads ITS1_R1/NICU_ITS1_R1only_sequences.qza \
  --o-classification ITS1_R1/NICU_ITS1_R1only_sklearn-UNITE-taxonomy.qza
qiime metadata tabulate \
  --m-input-file ITS1_R1/NICU_ITS1_R1only_sklearn-UNITE-taxonomy.qza \
  --o-visualization ITS1_R1/NICU_ITS1_R1only_sklearn-UNITE-taxonomy.qzv


###Use Miseq 300PE both reads and perform paired-end read joining

#bac16SV3V4
FAST/fast.py -generate_mapping -i bac16SV3V4/read1 -o bac16SV3V4/read1_map.txt
FAST/fast.py -generate_mapping -i bac16SV3V4/read2 -o bac16SV3V4/read2_map.txt
FAST/fast.py -add_labels -m bac16SV3V4/read1_map.txt -i bac16SV3V4/read1 -o bac16SV3V4/read1_labeled -t 2
FAST/fast.py -add_labels -m bac16SV3V4/read2_map.txt -i bac16SV3V4/read2 -o bac16SV3V4/read2_labeled -t 2
FAST/fast.py -merge_seqs -i bac16SV3V4/read1_labeled -o bac16SV3V4/read1.fastq
FAST/fast.py -merge_seqs -i bac16SV3V4/read2_labeled -o bac16SV3V4/read2.fastq
FAST/cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -o bac16SV3V4/read1.cut.fastq -p bac16SV3V4/read2.cut.fastq bac16SV3V4/read1.fastq bac16SV3V4/read2.fastq -m 50
rm -rf bac16SV3V4/read1_labeled/ bac16SV3V4/read2_labeled/
pear -f bac16SV3V4/read1.cut.fastq -r bac16SV3V4/read2.cut.fastq -o bac16SV3V4/merge.pear -k -j 2
mv bac16SV3V4/merge.pear.assembled.fastq bac16SV3V4/bac16SV3V4_merged.fastq
FAST/cutadapt -g CCTACGGGNGGCWGCAG...ATTAGAWACCCBNGTAGTCC -e 0.05 -o bac16SV3V4/bac16SV3V4_merged_trimmed.fastq bac16SV3V4/bac16SV3V4_merged.fastq --discard-untrimmed
rm -rf bac16SV3V4/read*fastq bac16SV3V4/merge.pear.* bac16SV3V4/bac16SV3V4_merged.fastq
FAST/vsearch --fastq_filter bac16SV3V4/bac16SV3V4_merged_trimmed.fastq --fastq_maxee 1 --fastaout bac16SV3V4/bac16SV3V4_merged_trimmed_filtered.fasta --fastq_minlen 250 --fasta_width 0
FAST/fast.py -stat_seqs -i bac16SV3V4/bac16SV3V4_merged_trimmed_filtered.fasta -o bac16SV3V4/bac16SV3V4_merged_trimmed_filtered.report.txt
FAST/fast.py -dereplicate -sizeout -i bac16SV3V4/bac16SV3V4_merged_trimmed_filtered.fasta -o bac16SV3V4/bac16SV3V4_merged_trimmed_filtered_derep -t 2
FAST/usearch -unoise3 bac16SV3V4/bac16SV3V4_merged_trimmed_filtered_derep.fasta -tabbedout bac16SV3V4/bac16SV3V4_unoise-denoise.txt -unoise_alpha 5
grep "\tchfilter\tzotu" bac16SV3V4/bac16SV3V4_unoise-denoise.txt | cut -f 1 | cut -d ";" -f 1 > bac16SV3V4/temp-list.txt
FAST/fast.py -filter_otu_map -i bac16SV3V4/bac16SV3V4_merged_trimmed_filtered_derep.txt -name_list bac16SV3V4/temp-list.txt -o bac16SV3V4/bac16SV3V4_merged_trimmed_filtered_derep_zotu.txt
FAST/fast.py -pick_seqs -sizeout -i bac16SV3V4/bac16SV3V4_merged_trimmed_filtered_derep.fasta -map bac16SV3V4/bac16SV3V4_merged_trimmed_filtered_derep_zotu.txt -o bac16SV3V4/bac16SV3V4_merged_trimmed_filtered_derep_zotu.fasta
FAST/fast.py -make_otu_table -qiime_map bac16SV3V4/bac16SV3V4_merged_trimmed_filtered_derep_zotu.txt -o bac16SV3V4/bac16SV3V4_otu-table.txt
FAST/vsearch --usearch_global bac16SV3V4/bac16SV3V4_merged_trimmed_filtered_derep_zotu.fasta -db ../rDNA_database/silva_nr_v132_train_set.fa --userout bac16SV3V4/bac16SV3V4_vsearch_taxa.txt --userfields query+target+ql+pairs+id --id 0.6
FAST/fast.py -filter_taxonomy -i bac16SV3V4/bac16SV3V4_vsearch_taxa.txt -op bac16SV3V4/temp.txt -match_length 0.9 -pident 97
perl -pe 's/;size=\d*\t/\t/' bac16SV3V4/temp.txt > bac16SV3V4/bac16SV3V4_vsearch_taxa.txt
FAST/fast.py -assign_taxonomy -keep_all -otu bac16SV3V4/bac16SV3V4_otu-table.txt -tax bac16SV3V4/bac16SV3V4_vsearch_taxa.txt -o bac16SV3V4/bac16SV3V4_otu-table_taxa.txt -scores
perl -pe 's/;size=\d*//' bac16SV3V4/bac16SV3V4_merged_trimmed_filtered_derep_zotu.fasta > bac16SV3V4/temp
mv bac16SV3V4/temp bac16SV3V4/bac16SV3V4_merged_trimmed_filtered_derep_zotu.fasta
rm -rf bac16SV3V4/temp*

#ITS1
FAST/fast.py -generate_mapping -i ITS1/read1 -o ITS1/read1_map.txt
FAST/fast.py -generate_mapping -i ITS1/read2 -o ITS1/read2_map.txt
FAST/fast.py -add_labels -m ITS1/read1_map.txt -i ITS1/read1 -o ITS1/read1_labeled -t 2
FAST/fast.py -add_labels -m ITS1/read2_map.txt -i ITS1/read2 -o ITS1/read2_labeled -t 2
FAST/fast.py -merge_seqs -i ITS1/read1_labeled -o ITS1/read1.fastq
FAST/fast.py -merge_seqs -i ITS1/read2_labeled -o ITS1/read2.fastq
FAST/cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -o ITS1/read1.cut.fastq -p ITS1/read2.cut.fastq ITS1/read1.fastq ITS1/read2.fastq -m 50
rm -rf ITS1/read1_labeled/ ITS1/read2_labeled/
pear -f ITS1/read1.cut.fastq -r ITS1/read2.cut.fastq -o ITS1/merge.pear -k -j 2
mv ITS1/merge.pear.assembled.fastq ITS1/ITS1_merged.fastq
FAST/cutadapt -g CTTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATT...AACTTTCAACAACGGATCTCTTGGYTCTSGCATCGATGAAGAACGCAGC -o ITS1/ITS1_merged_trimmed.fastq ITS1/ITS1_merged.fastq --discard-untrimmed
rm -rf ITS1/read*fastq ITS1/ITS1_merged.fastq ITS1/merge.pear.* 
FAST/vsearch --fastq_filter ITS1/ITS1_merged_trimmed.fastq --fastq_maxee 1 --fastaout ITS1/ITS1_merged_trimmed_filtered.fasta --fasta_width 0
FAST/fast.py -stat_seqs -i ITS1/ITS1_merged_trimmed_filtered.fasta -o ITS1/ITS1_merged_trimmed_filtered.report.txt
FAST/fast.py -dereplicate -sizeout -i ITS1/ITS1_merged_trimmed_filtered.fasta -o ITS1/ITS1_merged_trimmed_filtered_derep -t 2
FAST/usearch -unoise3 ITS1/ITS1_merged_trimmed_filtered_derep.fasta -tabbedout ITS1/ITS1_unoise-denoise.txt -unoise_alpha 5
grep "\tchfilter\tzotu" ITS1/ITS1_unoise-denoise.txt | cut -f 1 | cut -d ";" -f 1 > ITS1/temp-list.txt
FAST/fast.py -filter_otu_map -i ITS1/ITS1_merged_trimmed_filtered_derep.txt -name_list ITS1/temp-list.txt -o ITS1/ITS1_merged_trimmed_filtered_derep_zotu.txt
FAST/fast.py -pick_seqs -sizeout -i ITS1/ITS1_merged_trimmed_filtered_derep.fasta -map ITS1/ITS1_merged_trimmed_filtered_derep_zotu.txt -o ITS1/ITS1_merged_trimmed_filtered_derep_zotu.fasta
FAST/fast.py -make_otu_table -qiime_map ITS1/ITS1_merged_trimmed_filtered_derep_zotu.txt -o ITS1/ITS1_otu-table.txt
FAST/vsearch --usearch_global ITS1/ITS1_merged_trimmed_filtered_derep_zotu.fasta -db ../rDNA_database/sh_general_release_dynamic_01.12.2017.fasta --userout ITS1/ITS1_vsearch_taxa.txt --userfields query+target+ql+pairs+id --id 0.6
FAST/fast.py -filter_taxonomy -i ITS1/ITS1_vsearch_taxa.txt -op ITS1/temp.txt -match_length 0.9 -pident 97
perl -pe 's/;size=\d*\t/\t/' ITS1/temp.txt > ITS1/ITS1_vsearch_taxa.txt
FAST/fast.py -assign_taxonomy -keep_all -otu ITS1/ITS1_otu-table.txt -tax ITS1/ITS1_vsearch_taxa.txt -o ITS1/ITS1_otu-table_taxa.txt -scores
perl -pe 's/;size=\d*//' ITS1/ITS1_merged_trimmed_filtered_derep_zotu.fasta > ITS1/temp
mv ITS1/temp ITS1/ITS1_merged_trimmed_filtered_derep_zotu.fasta
rm -rf ITS1/temp*

#OTU assignment using qiime2
source activate qiime2-2018.2
qiime tools import \
  --input-path bac16SV3V4/bac16SV3V4_merged_trimmed_filtered_derep_zotu.fasta \
  --output-path bac16SV3V4/bac16SV3V4_sequences.qza \
  --type 'FeatureData[Sequence]'
qiime feature-classifier classify-sklearn \
  --i-classifier ../rDNA_database/silva-119-99-nb-classifier.qza \
  --i-reads bac16SV3V4/bac16SV3V4_sequences.qza \
  --o-classification bac16SV3V4/bac16SV3V4_sklearn-silva-taxonomy.qza
qiime metadata tabulate \
  --m-input-file bac16SV3V4/bac16SV3V4_sklearn-silva-taxonomy.qza \
  --o-visualization bac16SV3V4/bac16SV3V4_sklearn-silva-taxonomy.qzv

qiime tools import \
  --input-path ITS1/ITS1_merged_trimmed_filtered_derep_zotu.fasta \
  --output-path ITS1/ITS1_sequences.qza \
  --type 'FeatureData[Sequence]'
qiime feature-classifier classify-sklearn \
  --i-classifier ../rDNA_database/UNITE_s_01.12.2017_classifier.qza \
  --i-reads ITS1/ITS1_sequences.qza \
  --o-classification ITS1/ITS1_sklearn-silva-taxonomy.qza
qiime metadata tabulate \
  --m-input-file ITS1/ITS1_sklearn-silva-taxonomy.qza \
  --o-visualization ITS1/ITS1_sklearn-silva-taxonomy.qzv


