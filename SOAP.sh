#!/bin/bash

USAGE="Usage: bash $0 [options] <ref.fas> <fwd.fastq rvs.fastq>...

Options:
	-h          Prints help info
	-a          Attempt to retrieve abundance data [off]
	-b <int>    Kmer beginning value [35]
	-d <dir>    Directory where scripts exist [$(pwd)]
	-e <int>    Kmer ending value [95]
	-g <file>   File containing genera to check Blast output [None]
	-s <int>    Kmer stepping value [10]
	-t <int>    Number of threads [1]
"

clean_up() {
	#Deletes broken files if they were created and ends program
	printf "Creating $1" >&2
	rm -f $1
	shift
	for file in $@
	do
		printf ", $file" >&2
		rm -f $file
	done
	printf " failed\n" >&2
	exit 1
}

abundance=0
genera=""
kmerBegin=35
kmerEnd=95
kmerStep=10
scriptdir=$(pwd)
threads=1
seqs=0

while getopts ':hab:d:e:g:s:t:' flag
do
	case $flag in
		a) # Toggle abundance
			abundance=1
			;;
		h) # Help
			echo "$USAGE"
			exit 1
			;;
		b) # Kmer beginning value
			kmerBegin=$OPTARG
			;;
		d) # Change script directory
			scriptdir=$OPTARG
			;;
		e) # Kmer ending value
			kmerEnd=$OPTARG
			;;
		g) # File containing genera
			genera=$OPTARG
			;;
		s) # Kmer stepping value
			kmerStep=$OPTARG
			;;
		t) # Number of threads
			threads=$OPTARG
			;;
		\?)# Unexpected flag
			echo "Unexpected option -$OPTARG" >&2
			exit 1
			;;
		:) # Missing argument
			echo "Option -$OPTARG requires an argument" >&2
			exit 1
			;;
		*) # Unexpected error
		    echo "Unexpected error in getopts" >&2
			exit 1
			;;
	esac
done
shift $( expr $OPTIND - 1 )

if [ $# -le 1 ]
then
	echo $USAGE
	exit 1
fi

date +"Time started: %F %T"
echo "Using Reference Genome: $1"
ref=$1
shift
echo "Using $1 threads"

if [ $# -gt 0 ]
then
	#Create BLAST DB
	dbname="${ref%.*}_DB"
	dbdir=$(dirname $dbname)
	dbname=$(basename $dbname)
	echo $dbdir/$dbname
	if [ -f $dbdir/$dbname.nsq ]; then
		echo "BLAST database of the reference exists, skipping"
	else
		echo "Creating BLAST database"
		makeblastdb -in $ref -out $dbdir/$dbname -dbtype nucl || clean_up $dbdir/$dbname{.nhr,.nin,.nsq}
	fi
	
	#If using '..' or '.' to define relative location, then ref's location is augmented
	tmp=$( dirname $ref | cut -d / -f 1 )
	dirs=("." "..")
	dirs+=( $(ls -F | grep ".*/" | sed 's|/||') )
	for d in ${dirs[@]}
	do
		if [ "$tmp" = "$d" ]
		then
			ref="../../$ref"
			dbdir="../../$dbdir"
			genera="../../$genera"
			break
		fi
	done
	
	#Creating a new dir to store
	tmp=$(basename $ref)
	refname=${tmp%.*}	
	echo "Making ${refname}_SOAP directory"
	if [ ! -d $refname\_SOAP ]
	then
	  mkdir $refname\_SOAP
	fi
	cd $refname\_SOAP || { echo "Changing Directory failed"; exit 1; }
fi

while test ${#} -gt 0
do
	date +"Pair start: %T"
    seqname=$(basename $1)
	seqname=${seqname%.*}
	seqname=${seqname%_R*}
	echo "Using Sequence Name: $seqname"

	if [ ! -d "$seqname" ]
	then
		mkdir $seqname
	fi
	echo "Changing to directory $seqname"
	cd $seqname || { echo "Changing Directory failed"; exit 1; }
	
	fwd=$1
	echo "Forward read: $fwd"
	shift
	rvs=$1
	echo "Reverse read: $rvs"
	shift

	#Run de novo assembler with different kmer values
	if [ -f $seqname\_SOAP.fasta ]; then
		echo "SOAP alread ran, skipping"
	else
		for ((i=$kmerBegin; i <= $kmerEnd; i += $kmerStep)); do
			if [ -d SOAP_Run_$i ]; then
				echo "SOAP with kmer $i already ran, skipping"
			elif [ -f SOAP_Run_$i.tar.gz ]; then
				echo "Unzipping SOAP_Run_$i.tar.gz"
				tar xfz SOAP_Run_$i.tar.gz
			else
				mkdir SOAP_Run_$i
				date +"SOAP_Run_$i start: %T"
				echo "Creating SOAP_Run_$i/$seqname.config"
				cat <<EOT > SOAP_Run_$i/$seqname.config
#maximal read length
max_rd_len=150
[LIB]
#maximal read length in this lib
rd_len_cutoff=150
#average insert size
avg_ins=275
#if sequence needs to be reversed
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
#fastq file for read 1
q1=$fwd
#fastq file for read 2 always follows fastq file for read 1
q2=$rvs
EOT
				echo "Running SOAP with kmer: $i"
				~/SOAPdenovo-Trans/SOAPdenovo-Trans-127mer all -s SOAP_Run_$i/$seqname.config -K $i -p $threads -o SOAP_Run_$i/$seqname || { echo "Running SOAP failed"; exit 1; }
				echo "Done"
				date +"SOAP_Run_$i stop: %T"
			fi
			python $scriptdir/PrefixFasta.py -p "SOAP_Run_$i|" -r " " "_" -o SOAP_Run_$i.fasta SOAP_Run_$i/$seqname.contig SOAP_Run_$i/$seqname.scafSeq || clean_up SOAP_Run_$i.fasta $seqname\_SOAP.fasta
			cat SOAP_Run_$i.fasta >> $seqname\_SOAP.fasta
			rm -f SOAP_Run_$i.fasta
			if [ ! -f SOAP_Run_$i.tar.gz ]; then
				tar cfz SOAP_Run_$i.tar.gz SOAP_Run_$i
			fi
			rm -R SOAP_Run_$i
		done
	fi

	#Remove redundancies between multiple kmer runs
	if [ -f $seqname\_SOAP_NR.fasta ]; then
		echo "Redundancies removed, skipping"
	else
		echo "Removing redundancies"
		cd-hit-est -c 0.9 -n 8 -r 1 -T $threads -M 0 -i $seqname\_SOAP.fasta -o $seqname\_SOAP_NR.fasta || clean_up $seqname\_SOAP_NR.fasta $seqname\_SOAP_NR.fasta.clstr
		rm -f $seqname\_SOAP_NR.fasta.clstr
	fi
	
	#Collect N50, L50 and other related statistics
	if [ -f NandLStats_$seqname\_SOAP_NR.txt ]; then
		echo "N and L stats calculated, skipping"
	else
		echo "Calculating N and L stats"
		perl $scriptdir/NandLStats.pl $seqname\_SOAP_NR.fasta > NandLStats_$seqname\_SOAP_NR.txt || clean_up NandLStats_$seqname\_SOAP_NR.txt
	fi
	
	#Blast contigs
	if [ -f blasted_$seqname\_SOAP_$dbname.xml ]; then
		echo "SOAP output blasted, skipping"
	else
		echo "Blasting SOAP output"
		blastn -db $dbdir/$dbname -query $seqname\_SOAP_NR.fasta -outfmt 5 -task megablast -max_target_seqs 5 -word_size 7 -evalue 0.00001 -num_threads $threads -out blasted_$seqname\_SOAP_$dbname.xml || clean_up blasted_$seqname\_SOAP_$dbname.xml
	fi

	#Create transcriptome via scaffolding
	if [ -d SOAP_Spliced ]; then
		echo "Transcriptome assembled, skipping"
	else
		echo "Assembling transcriptome"
		python $scriptdir/AnnotateContigsV1.?.py -C SOAP_Contigs -s SOAP_Spliced -r -c 0 -t 0.7 -R $ref -x blasted_$seqname\_SOAP_$dbname.xml -m cds || { echo "AnnotateContigsV1.? failed"; exit 1; }
	fi

	#Condense transcriptome into one file
	if [ -f $seqname\_DeNovo_All_SOAP.fasta ]; then
		echo "Transcriptome condensed into one file, skipping"
	else
		echo "Condensing transcriptome into one file"
		python $scriptdir/ParseMerged.py $seqname\_DeNovo_All_SOAP.fasta $(ls SOAP_Spliced/*.fasta) || clean_up $seqname\_DeNovo_All_SOAP.fasta
	fi

	if [ $abundance -gt 0 ]; then
		#Prepare for creating sam file
		if [ -d RSEM ]; then
			echo "RSEM dir created, skipping"
		else
			echo "Creating RSEM dependencies"
			mkdir RSEM
			rsem-prepare-reference --bowtie $seqname\_SOAP_NR.fasta RSEM/$seqname\_SOAP_NR || { echo "Preparing RSEM dependencies failed"; rm -R RSEM; exit 1; }
		fi
		
		#Create sam file
		if [ -f RSEM/$seqname\_SOAP_NR.sam ]; then
			echo "SAM exists, skipping"
		else
			echo "Creating SAM"
			bowtie -p $threads --all -S RSEM/$seqname\_SOAP_NR -1 $fwd -2 $rvs > RSEM/$seqname\_SOAP_NR.sam || clean_up RSEM/$seqname\_SOAP_NR.sam
		fi
		
		#Create abundance data
		if [ -f $seqname\_SOAP_NR.genes.results ]; then
			echo "Abundance data alread created, skipping"
		else
			echo "Creating abundance data"
			rsem-calculate-expression -p $threads --paired-end $fwd $rvs RSEM/$seqname\_SOAP_NR $seqname\_SOAP_NR || clean_up $seqname\_SOAP_NR.{genes.results,isoforms.results,transcripts.bam}
			rm -f $seqname\_SOAP_NR.transcripts.bam
			rm -R $seqname\_SOAP_NR.stat
		fi
		
		#Condense isoform abundance data
		if [ -f annotated_isoforms_$seqname.xlsx ]
		then
			echo "Annotation by isoforms exists, skipping"
		else
			echo "Annotating by isoforms"
			python $scriptdir/TranscriptAnnotate.py blasted_$seqname\_SOAP_$dbname.xml $seqname\_SOAP_NR.isoforms.results annotated_isoforms_$seqname.xlsx || clean_up annotated_isoforms_$seqname.xlsx
		fi
		
		#Condense gene abundance data
		if [ -f annotated_genes_$seqname.xlsx ]
		then
			echo "Annotation by isoforms exists, skipping"
		else
			echo "Annotating by isoforms"
			python $scriptdir/TranscriptAnnotateByGenes.py blasted_$seqname\_SOAP_$dbname.xml $seqname\_SOAP_NR.genes.results annotated_genes_$seqname.xlsx || clean_up annotated_genes_$seqname.xlsx
		fi
	fi
		
	date +"Pair end: %T"
	echo "_________________________________________________"
	echo "FINISHED. Moving on to next pair"
	echo "_________________________________________________"
	cd ..  || { echo "Changing Directory failed"; exit 1; }
done

#Condense Summary_Results.csv from each taxon
if [ -f Summary.csv ]; then
	echo "Breadth of coverage results condensed, skipping"
else
	echo "Condensing breadth of coverage results"
	python $scriptdir/DeNovoSummary.py . Summary.csv
fi

echo "All Done here!"
date +"Time ended: %F %T"
