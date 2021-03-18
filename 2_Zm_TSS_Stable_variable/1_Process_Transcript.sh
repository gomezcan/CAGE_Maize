#/bin/bash

# 1. Get transcripts
grep "ID=transcript:" Zea_mays.B73_RefGen_v4.46.gff3 | cut -f9 | sed $'s/;/\t/g' | sed 's/ID=transcript://g' | sed 's/Parent=gene://g' \
	| sed 's/biotype=//g' | sed 's/transcript_id=//g' | cut -f2-4 | sort -u -o Zm_Gene_Transcript.txt ;

# 2. Summary transcript by gene
cut -f1 Zm_Gene_Transcript.txt | sort | uniq -c | awk '{ print $2,$1}' > Zm_Freq.Transcript.by.Gene.txt;

# 3. Get annotated TSS

awk -v OFS='\t' '{ if($3=="gene") print $0}'<  Zea_mays.B73_RefGen_v4.46.gff3 | sed $'s/;/\t/g' | cut -f4,5,7,9 | \
	awk -v OFS='\t' '{if($3=="+") print $4,$1,$3; else print $4,$2,$3}' | sed 's/ID=gene://g' > Zm_Annotated_TSS.txt  

# 4. Start position of first exon
### if + get the first coordenate, and if -  get the second coordenate

grep '.exon2' Zea_mays.B73_RefGen_v4.46.gff3 | sed $'s/;/\t/g' | sed 's/Parent=transcript://g' | cut -f4,5,7,9 | 
	awk -v OFS='\t' '{if($3=="+") print $4,$1,$3; else print $4,$2,$3}' > Zm_Exon2.Start_by_Transcript.txt

