#!/bin/sh

OPTIONS=$(getopt -o f:r:n:c:h --long pbdagcon,dom -n 'parse-options' -- "$@")

if [ $? -ne 0 ]; then
  echo "getopt error"
  exit 1
fi

echo "$OPTIONS"
eval set -- $OPTIONS

TOOL=0
NBSIMU=1
COVERAGE=100
while true; do
  case "$1" in
    
    -f) FILE="$2" ; shift ;;
    -r) REF="$2"; shift ;;
    -n) NBSIMU="$2";  shift ;;
    -c) COVERAGE="$2";  shift ;;
    --pbdagcon) TOOL=1 ;;
    --dom) TOOL=2;;
    -h) HELP=1 ;;
    --)        shift ; break ;;
    *)         echo "unknown option: $1" ; exit 1 ;;
  esac
  shift
done

if [ $# -ne 0 ]; then
  echo "unknown option(s): $@"
  exit 1
fi

echo "read file: $FILE"
echo "reference sequence: $REF"


REFLENGTH=$(tail -1 "$REF" | wc -c)
R --no-save <<RSCRIPT
df <- data.frame("consensus_mean_size", "simulated_reads_mean_size", "mean_backbone_similarity", "mean_consensus_similarity")
write.table(df, "stats.txt", row.names=F, col.names=F, quote=F)
RSCRIPT

nbsimu=0
#"${readcount} ${sizeB} ${sizeC} ${REFLENGTH} ${covered} ${simil_backbone} ${simil}"
echo "1 2 3 4 5 6 7 8 9" > results_similarity.txt
#~ echo "read_number_for_backbone size_backbone size_consensus size_ref percentage_of_ref_covered  similarity_backbone similarity_consensus" > results_similarity.txt

while (($nbsimu < $NBSIMU))
do
	echo "***** Simulation $nbsimu *****"
	echo "simulating data in ${FILE}_0001.fasta ..."
	/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/Linux-amd64/bin/pbsim $REF  --prefix $FILE --depth $COVERAGE  --length-mean $REFLENGTH --length-min 1000 --model_qc ~/bin/pbsim-1.0.3_modified/data/model_qc_ccs 2> /tmp/pbsim.log
	awk '{if ( NR%2==0 ){print length ($0)}}' ${FILE}_0001.fasta >> /tmp/reads_sizes.list
	NUMOFLINES=$(wc -l < "${FILE}_0001.fasta")
	i=1
	readcount=1
	#~ echo "read_number_for_backbone size_consensus size_ref percentage_of_ref_covered similarity_of_covered_ref_and_consensus" > results_similarity.txt

	if [ $TOOL = 1 ]; then
		echo "computing consensus, using PBDAGCON ..."
	else if [ $TOOL = 2 ]; then
			echo "computing consensus, using dom script ..."
		else
			echo "computing consensus, using SparC ..."
		fi
	fi

	while (($i < $NUMOFLINES))
	do
		j=$((i+1))
		sed -n "${i},${j}p" ${FILE}_0001.fasta > backbone_${i}.fasta
		sizeB=$(sed -n "2p" backbone_${i}.fasta | wc -c) #size of the backbone
		awk '{if ( NR%4==1 || NR%4==2){print $0}}' ${FILE}_0001.fastq  | sed 's/@/>/g' > ${FILE}_0001.fasta
		~/bin/blasr -nproc 2 backbone_${i}.fasta ${REF} -bestn 1 -m 5 -minMatch 15 -out blasr_bb.mapped.m5
		matchesbb=$(awk '{print $12}' blasr_bb.mapped.m5)
		nbIns=$(awk '{print $14}' blasr_bb.mapped.m5)
		nbDel=$(awk '{print $15}' blasr_bb.mapped.m5)
		nbIndel=$((nbIns+nbDel))
		#qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand score numMatch numMismatch numIns numDel mapQV qAlignedSeq matchPattern tAlignedSeq
		qstart=$(awk '{print $3}' blasr_bb.mapped.m5)
		qend=$(awk '{print $4}' blasr_bb.mapped.m5)
		tstart=$(awk '{print $8}' blasr_bb.mapped.m5)
		tend=$(awk '{print $9}' blasr_bb.mapped.m5)
		shift=$(( (qstart-tstart)*(qstart-tstart)+(qend-tend)*(qend-tend) ))
		#~ echo $shift
		if [ "$matchesbb" = "" ]; then
			matchesbb=0
		fi
		simil_backbone=$(awk "BEGIN {printf \"%.2f\",${matchesbb}/${sizeB}}")
		rm blasr_bb.mapped.m5
		~/bin/blasr -nproc 2 ${FILE}_0001.fasta backbone_${i}.fasta -bestn 1 -m 5 -minMatch 15 -out blasr_output.mapped.m5
		if [ $TOOL = 1 ]; then
			~/bin/pbdagcon/src/cpp/pbdagcon blasr_output.mapped.m5 2> /tmp/consensus.log > Consensus.consensus.fasta
		else if [ $TOOL = 2 ]; then
				tail -1 backbone_${i}.fasta > tmp_reads.fa
				awk '{if ( NR%2==0 ){print $0}}' ${FILE}_0001.fasta >> tmp_reads.fa
				./correctLR.py tmp_reads.fa > tmp_consensus.fa
				echo ">c" > Consensus.consensus.fasta
				tail -1 tmp_consensus.fa >> Consensus.consensus.fasta
			else 
				~/bin/Sparc_Linux/Sparc m blasr_output.mapped.m5  b backbone_${i}.fasta  k 3 c 3 g 3 > /tmp/consensus.log
			fi
		fi
		~/bin/blat $REF Consensus.consensus.fasta blat.psl > /tmp/blat.log 2> /dev/null
		sizeC=$(sed -n "2p" /tmp/blat.log | awk '{print $2}') #size of the consensus
		if [ "$sizeC" = "" ]; then
			sizeC=0
			covered=0
			match=0
			simil=0
		else
			covered=$(awk "BEGIN {printf \"%.2f\", ${sizeC}/${REFLENGTH}}")
			match=$(sed -n "6p" blat.psl | awk '{print $1}')
			if [ "$match" = "" ]; then
				match=0
			fi
			simil=$(awk "BEGIN {printf \"%.2f\",${match}/${sizeC}}")
		fi
		echo "${readcount} ${sizeB} ${sizeC} ${REFLENGTH} ${covered} ${simil_backbone} ${simil} ${nbIndel} ${shift}" >> results_similarity.txt
		((i++))
		((i++))
		((readcount++))
		 
	done
	echo "writing results in results_similarity.txt"

	rm backbone_*
	rm blasr_output.mapped.m5
	if [ $TOOL = 2 ]; then
		rm tmp_reads.fa
		rm tmp_consensus.fa
	else
		rm  /tmp/consensus.log
	fi

	((++nbsimu))
done

R --no-save <<RSCRIPT
pdf("distribution_of_consensus_size.pdf")
res <- read.table("results_similarity.txt", h=T, sep=" ")
sizes <- read.table("/tmp/reads_sizes.list")
hist(res[,3], xlab="size of consensus", ,main="distribution of consensus sizes")
df <- data.frame(mean(res[,3]), mean(sizes[,1]), mean(res[,7]))
names(df) <- c("consensus_mean_size", "simulated_reads_mean_size", "mean_similarity")
write.table(df, "stats.txt", row.names=F, col.names=F, quote=F, append=T)
dev.off()
pdf("consensus_and_similarity_of_consensus_in_function_of_backbone.pdf")
par(mfrow = c(3, 3))
plot(res[,2]~res[,3], xlab="consensus_size", ylab="backbone_size")
plot(res[,6]~res[,7], xlab="consensus_similarity", ylab="backbone_similarity")
plot(res[,2]~res[,7], xlab="consensus_similarity", ylab="backbone_size")
plot(res[,6]~res[,3], xlab="consensus_size", ylab="backbone_similarity")
plot(res[,8]~res[,3], xlab="consensus_size", ylab="backbone_indels")
plot(res[,8]~res[,7], xlab="consensus_similarity", ylab="backbone_indels")
plot(res[,9]~res[,7], xlab="consensus_similarity", ylab="backbone_shift")
plot(res[,9]~res[,3], xlab="consensus_size", ylab="backbone_shift")
dev.off()
RSCRIPT

rm /tmp/reads_sizes.list

#qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand score numMatch numMismatch numIns numDel mapQV qAlignedSeq matchPattern tAlignedSeq
