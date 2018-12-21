#!/bin/bash

#clip and filter reads
clip.sh

#merge paired reads
usearch -fastq_mergepairs clip/*_R1.fastq -fastqout merged.fastq -relabel @ -fastq_minovlen 10 -minhsp 8 -fastq_maxdiffpct 20

#Filter merged reads below 300 bp
bysize.py merged.fastq sized.fastq fastq 300 20000
#convert to fasta
convert.py sized.fastq fastq fasta

#UPARSE
usearch -derep_fulllength sized.fasta -sizeout -fastaout uniques.fa
usearch -cluster_otus uniques.fa -minsize 2 -otus otus.fa -relabel Otu
usearch -utax otus.fa -db ~/db/utax/16s.udb -strand both -strand both -fastaout otus_tax.fa -utax_cutoff 0.9 -log utax.log
usearch -usearch_global otus.fa -db ~/db/utax/16s.udb -strand both -id 0.97 -alnout otus_ref.aln -userout otus_ref.user -userfields query+target+id
usearch -usearch_global sized.fasta -db otus_tax.fa -strand both -id 0.97 -log make_otutab.log -otutabout otutab.txt -biomout otutab.biom

#Convert biom to text for viewing
biom convert -i otutab.biom -o otutab.txt --to-tsv

#rarefaction
alpha_rarefaction.py -i otutab.biom -m ../47.1/mapping.txt -o rarefaction1 -n 25 -f -p ../47.1/qiime_parameters1.txt

#0.5% abundance filter, remove singleton samples' otus:
filter_otus_from_otu_table.py -i otutab.biom -o otufilt.biom --min_count_fraction 0.001 -s 2
normalize_table.py -i otufilt.biom -o otunorm.biom

biom convert -i otufilt.biom -o otufilt.txt --to-tsv --header-key taxonomy && sed -i 's/; //g' otufilt.txt && sed -i 's/"//g' otufilt.txt
biom convert -i otunorm.biom -o otunorm.txt --to-tsv --header-key taxonomy && sed -i 's/; //g' otunorm.txt && sed -i 's/"//g' otunorm.txt

#get list of filtered otus and extract to new fasta
cut -f1 otufilt.txt | sed -n '1!p' > otu.list
get-seqs-from-file.py otu.list otus.fa otufilt.fa fasta fasta

#draw tree
muscle -in otufilt.fa -out otusfilt.aln
FastTree otusfilt.aln > otusfilt.tree

#calculate unifrac matrix
#sort by view order!
sort_otu_table.py -i otunorm.biom -o otusort.biom -m mapping.txt -s View_order
beta_diversity.py -i otusort.biom -m weighted_unifrac -o ./ -t otusfilt.tree

#rarefaction and alpha on filtered table
alpha_rarefaction.py -i otufilt.biom -m mapping.txt -o rarefaction2 -n 25 -f -p qiime_parameters1.txt

#compare alpha diversity n.b. set depth to lowest abundance sample
compare_alpha_diversity.py -i rarefaction2/alpha_div_collated/simpson.txt -m mapping.txt -c Diet,Drug,Description -o compare_simpson -t nonparametric -d 9178
compare_alpha_diversity.py -i rarefaction2/alpha_div_collated/chao1.txt -m mapping.txt -c Diet,Drug,Description -o compare_chao1 -t nonparametric -d 9178
compare_alpha_diversity.py -i rarefaction2/alpha_div_collated/shannon.txt -m mapping.txt -c Diet,Drug,Description -o compare_shannon -t nonparametric -d 9178
compare_alpha_diversity.py -i rarefaction2/alpha_div_collated/observed_species.txt -m mapping.txt -c Diet,Drug,Description -o compare_obs -t nonparametric -d 9178

#PCO
principal_coordinates.py -i weighted_unifrac_otusort.txt -o pco.out

#test differnce between treatment groups
group_significance.py -i otunorm.biom -m anova-map.txt -o sd-mino_anova.txt -c sd-mino-veh -s ANOVA 
group_significance.py -i otunorm.biom -m anova-map.txt -o hf-mino_anova.txt -c hf-mino-veh -s ANOVA 
group_significance.py -i otunorm.biom -m anova-map.txt -o sd-lnac_anova.txt -c sd-lnac-veh -s ANOVA 
group_significance.py -i otunorm.biom -m anova-map.txt -o hf-lnac_anova.txt -c hf-lnac-veh -s ANOVA 
group_significance.py -i otunorm.biom -m anova-map.txt -o sd-hnac_anova.txt -c sd-hnac-veh -s ANOVA 
group_significance.py -i otunorm.biom -m anova-map.txt -o hf-hnac_anova.txt -c hf-hnac-veh -s ANOVA 

#produce PCO plot of CSS normalised data
normalize_table.py -i otusort.biom -o css_norm.biom -s 
beta_diversity.py -i css_norm.biom -m weighted_unifrac -o ./ -t otusfilt.tree
principal_coordinates.py -i weighted_unifrac_css_norm.txt -o pco.out
