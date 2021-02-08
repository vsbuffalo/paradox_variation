all: dmel-all-r5.33.gff.gz dmel-all-chromosome-r5.33.fasta.gz dmel-5.33-seqlens.txt

dmel-all-r5.33.gff.gz: 
	curl -o  dmel-all-r5.33.gff.gz ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.33_FB2011_01/gff/dmel-all-r5.33.gff.gz

dmel-all-chromosome-r5.33.fasta.gz: 
	curl -o dmel-all-chromosome-r5.33.fasta.gz ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.33_FB2011_01//fasta/dmel-all-chromosome-r5.33.fasta.gz

dmel-5.33-seqlens.txt: dmel-all-chromosome-r5.33.fasta.gz
	bioawk -cfastx '{print $$name, length($$seq)}' dmel-all-chromosome-r5.33.fasta.gz > $@
