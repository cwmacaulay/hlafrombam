#!/bin/bash
# readphasing.sh
# Charles Macaulay
# 2022-03-18



# NOTE: BAM FILE NEEDS TO BE INDEXED. 
#BAMFILE=$1
#CONSTANT=$2
#REFGENOME=$3

# get unique HLAs for this bam file. 
#hlas=$( /N/u/cmacaulay/tools/bin/samtools view $BAMFILE | cut -f 3 | grep -v "*" | uniq )

# grab each one and pull it out of the bam file
#for hla in $( echo $hlas); do 
#readcenters=$(/N/u/cmacaulay/tools/bin/samtools view $BAMFILE "$hla" | awk '{ print $4 + (length($10)/2)}' | sort -n)
#num=$( /N/u/cmacaulay/tools/bin/samtools view $BAMFILE "$hla" | wc -l)
#lastcenter=$( echo "${readcenters%% *}" | head -1 );
#lastsum=0.0
#if (( num > 1 )); then
#    for center in $( echo $readcenters); do
#	dis=$( echo 'scale=4; '$center'-'$lastcenter'' | bc )
#	frac=$( echo 'scale=4; '$dis' / '$CONSTANT'' | bc )
#	add=$( echo 'scale=4; '$frac'*'$frac'*'$frac'' | bc )
#	sum=$( echo 'scale=4; '$lastsum'+'$add'' | bc )
#	lastsum=$sum
#	lastcenter=$center
#    done
#fi
#echo $hla $lastsum >> testout.txt
#done

# now need to make a reference genome using candidate alleles that have a
# score less than 125. 
#awk -F" " '$2 < 125' testout.txt | awk -F" " '$2 > 0.0' | awk -F " " '{ print $1 }' > hlainclude.txt
#cat $REFGENOME | sed -z 's/\n/bbb/g' | sed -z 's/>/\naaa/g' | grep -f /N/u/cmacaulay/20220125haplotyper/testreadphasing/hlainclude.txt | sed -z 's/\naaa/>/g' | sed -z 's/bbb/\n/g' | sed -z 's/aaa/>/g' > newgenome.fasta


#cat newgenome.fasta | sed -z 's/bp\n/bp\t/g' | sed -z 's/\n//g' | sed -z 's/>/\n>/g' > newgenomesearchable.txt

# to find and remove totally duplicate alleles. 
# get the duplicate allele sequences. 
cat "newgenomesearchable_10.txt" | awk -F " " '{ print $5 }' | sort | uniq -d > dups.txt
for i in $( cat dups.txt ); do
    id=$( cat "newgenomesearchable_10.txt" | grep $i | awk -F " " '{ print $1" "$2}' | head -1 | awk -F ":" '{ print $1":"$2":"$3 }' )
    cat newgenomesearchable_10.txt | grep -v $i > newgenomesearchable_10.nodup.tmp.txt
    echo $id $(echo $i | wc -c) bp $i >> newgenomesearchable_10.nodup.tmp.txt
    mv newgenomesearchable_10.nodup.tmp.txt newgenomesearchable_10.txt
done
cp newgenomesearchable_10.txt newgenomesearchableNODUPs.txt

# to find and remove 150bp repeated regions on shoulders of alleles. 
input="newgenomesearchable_10.txt"
while IFS= read -r line; do
    allele=$( echo $line | awk -F " " '{print $5}')
    id=$( echo $line | awk -F " " '{print $1}')
   
    # need to separate the genome file so that we aren't searching within a given allele for duplicated
    # sequences. 
    cat newgenomesearchable_10.txt | grep -v $id > genomesplit.txt
    cat genomesplit.txt >> outputmonitoring.txt
    declare -i WINSIZE=150
    declare -i SLIDING=1
    NELEM=$(echo $allele | wc -c )
    NSLIDES=$(((NELEM - WINSIZE) / SLIDING))
    for ((i=0; i<=NSLIDES;i++)); do
        START=$((1+$i*SLIDING))
        END=$((START+(WINSIZE)))
	search=$( echo $allele | awk '{ print substr( $0, '$START', $START+150 ) }' )
	cat genomesplit.txt | sed -s 's/ '$search'/ /g' > genomesplitresult.txt
	mv genomesplitresult.txt genomesplit.txt
	cat genomesplit.txt | sed -z 's/'$search'\n/\n/g' > genomesplitresult.txt
	mv genomesplitresult.txt genomesplit.txt
    done
    echo $line >> genomesplit.txt
    mv genomesplit.txt newgenomesearchable_10.txt
done < "$input"


# now need to get the positions to rule in:
input="newgenomesearchableNODUPs.txt"
while IFS= read -r line; do
    origallele=$( echo $line | awk -F " " '{print $5}')
    id=$( echo $line | awk -F " " '{print $1}')
    shortenedallele=$( cat newgenomesearchable_10.txt | grep $id | awk -F " " '{print $5}' )
    echo "~~~~~~~~~~~~~~"
    echo $origallele
    echo ">>>>>>>>>>>>>>"
    echo $shortenedallele
    echo "~~~~~~~~~~~~~~"
    rest=${origallele#*$shortenedallele}
    echo $(( ${#origallele} - ${#rest} - ${#searchstring} )) $( echo $shortenedallele | wc -c)
    echo "      "
done < "$input"
# NEED TO HANDLE HOLES IN THE MIDDLE!
