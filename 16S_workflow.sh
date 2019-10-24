#!/bin/bash
mkdir $2
mkdir $2/trimmed_output


#faire une boucle pour faire chaque paire ulterieurement
fastqc $1/*.fastq.gz -o $2/

gunzip $1/*.gz


#Trim les fastq avec AlienTrimmer, ranger les outputs dans un dossier output/trimmed_output
liste_fichiers=`ls $1/*_R1.fastq`

for fichier in $liste_fichiers
do
    echo $fichier
    R1="$fichier"
    R2=$(echo $fichier | sed s:R1:R2:g)
    echo $R2
    java -jar ./soft/AlienTrimmer.jar -if $R1 -ir $R2 -c ./databases/contaminants.fasta -q 20 -of ./$2/trimmed_output/$(basename $R1) -or ./$2/trimmed_output/$(basename $R2)
done

mkdir $2/merged
list_trimmed=`ls $2/trimmed_output/*_R1.fastq`
for fichier in $list_trimmed
do
    echo $fichier
    R1="$fichier"
    R2=$(echo $fichier | sed s:R1:R2:g)
    echo $R2
    name=$(echo $(basename $fichier) | cut -d. -f1)
    output=$name".merged.fastq"
    vsearch --fastq_mergepairs $R1 --reverse $R2 --fastaout $2/merged/$(basename $fichier).fasta --label_suffix ";sample=$name;"
done

#Merge everything into a single file
cat $2/merged/*.fasta | sed -e 's/ //g' > $2/amplicon.fasta

#Dereplication
mkdir $2/dereplicated
vsearch --derep_fulllength $2/amplicon.fasta --sizeout --output $2/dereplicated/derep.fasta

#Deleting singletons
vsearch --fastx_filter $2/dereplicated/derep.fasta --minsize 10 --fastaout $2/dereplicated/singsuppr.fasta

#Deleting chimeras
vsearch --uchime_denovo $2/dereplicated/singsuppr.fasta --chimeras $2/dereplicated/chimeras.fasta --nonchimeras $2/dereplicated/nonchim.fasta

#Clustering
vsearch --cluster_size $2/dereplicated/nonchim.fasta --id 0.97 --centroids $2/dereplicated/OTU.fasta --relabel OTU_


vsearch --usearch_global $2/amplicon.fasta --db $2/dereplicated/OTU.fasta --id 0.97 --otutabout $2/dereplicated/abundance_table.txt --sizeout

vsearch --usearch_global $2/amplicon.fasta --db databases/mock_16S_18S.fasta --id 0.9 --top_hits_only --userfields query+target --userout $2/dereplicated/vs_16S.txt

