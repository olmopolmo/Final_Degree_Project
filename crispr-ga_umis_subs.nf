#!/usr/bin/env nextflow

outDir = params.outDir

process mergeReads {
        
	publishDir "$outDir/Data", mode: 'copy'

	input:
	set val(sampleID), file(reads) from Channel.fromFilePairs( params.input, size : params.inSize )

	output:
	tuple sampleID, file("*assembled_t.fastq") into preprocessedReads
	tuple sampleID, file("*assembled_t.fastq") into noUMIpreprocessedReads
	tuple sampleID, file("read_0*.fastq"), file("*assembled.fastq"), file("*assembled_t.fastq"), file("*_adapters.out") into summaryRawMerge

	script:
	if ( params.r2file )
	"""
	gzip -dfc ${reads[0]} > read_0.fastq
	gzip -dfc ${reads[1]} > read_1.fastq
	pear -f read_0.fastq -r read_1.fastq -o $sampleID
	
	### Trimming adapters and quality filter
	fastqc ${sampleID}.assembled.fastq --extract
	test_lines=`grep -A 4 "Overrepresented sequences" ${sampleID}.assembled_fastqc/fastqc_data.txt | grep -v "No Hit" | head -3 | tail -1 | awk '{print NF}'`
	if [ \$test_lines -lt 6 ]; then 
		echo "Reads with adapters: 0 (0%)" > ${sampleID}_adapters.out
		fastq_quality_trimmer -Q33 -t 20 -l 80 -i ${sampleID}.assembled.fastq -o ${sampleID}.assembled_t.fastq; #trimming 3', 5'
	else 
		adapter_seq=`grep -A 4 "Overrepresented sequences" ${sampleID}.assembled_fastqc/fastqc_data.txt | grep -v "No Hit" | head -3 | tail -1 | awk -F "\t" '{ print \$1 }'`;
		cutadapt -g \$adapter_seq -o ${sampleID}.assembled_adapter-trim.fastq ${sampleID}.assembled.fastq > ${sampleID}_adapters.out
		fastq_quality_trimmer -Q33 -t 20 -l 80 -i ${sampleID}.assembled.fastq -o ${sampleID}.assembled_t.fastq; #trimming 3', 5' 
	fi
	
	"""
	else if ( params.simGE )
	"""
	touch read_0.dummy.fastq
	touch dummy.assembled.fastq
	touch dummy.assembled_t.fastq
	echo "Reads with adapters: 0 (0%)" > ${sampleID}_adapters.out
	"""
	else
	"""
	gzip -dfc ${reads[0]} > read_0.fastq
	cp read_0.fastq ${sampleID}.assembled.fastq
	
	### Trimming adapters and quality filter
	fastqc ${sampleID}.assembled.fastq --extract
	test_lines=`grep -A 4 "Overrepresented sequences" ${sampleID}.assembled_fastqc/fastqc_data.txt | grep -v "No Hit" | head -3 | tail -1 | awk '{print NF}'`
	if [ \$test_lines -lt 6 ]; then 
		echo "Reads with adapters: 0 (0%)" > ${sampleID}_adapters.out
		fastq_quality_trimmer -Q33 -t 20 -l 80 -i ${sampleID}.assembled.fastq -o ${sampleID}.assembled_t.fastq; #trimming 3', 5'
	else 
		adapter_seq=`grep -A 4 "Overrepresented sequences" ${sampleID}.assembled_fastqc/fastqc_data.txt | grep -v "No Hit" | head -3 | tail -1 | awk -F "\t" '{ print \$1 }'`;
		cutadapt -g \$adapter_seq -o ${sampleID}.assembled_adapter-trim.fastq ${sampleID}.assembled.fastq > ${sampleID}_adapters.out
		fastq_quality_trimmer -Q33 -t 20 -l 80 -i ${sampleID}.assembled.fastq -o ${sampleID}.assembled_t.fastq; #trimming 3', 5' 
	fi
	
	"""

}

process summary_merge{
	
	input:
	set sampleID, file(raw_reads), file(assembled_reads), file(trimmed_reads), file(trimmed_adapters) from summaryRawMerge
	
	output:
	tuple sampleID, file("*_summary.csv") into summaryReport
	
	script:
	"""
	#### Counts
	pre_raw_count=`cat ${raw_reads} | wc -l`;
	raw_count=`echo \$(( \$pre_raw_count / 4 ))`;
	pre_merged_count=`cat ${assembled_reads} | wc -l`;
	merged_count=`echo \$(( \$pre_merged_count / 4 ))`;
	adapters=`grep "Reads with adapters" ${trimmed_adapters} | awk -F " " '{ print \$(NF-1)" "\$NF }'`;
	pre_filt=`cat ${trimmed_reads} | wc -l`;
	filt=`echo \$(( \$pre_filt / 4 ))`;
	#### Percentages
	merged_perc=`echo "scale=1;(\$merged_count/\$raw_count)*100" | bc`
	filt_perc=`echo "scale=1;(\$filt/\$merged_count)*100" | bc`
	#### Table
	echo "class, count" > ${sampleID}_summary.csv
	echo "raw-reads," \$raw_count "(100.0%)" >> ${sampleID}_summary.csv
	echo "merged-reads," \$merged_count "("\$merged_perc"%)" >> ${sampleID}_summary.csv
	echo "reads-with-adapters," \${adapters//,} >> ${sampleID}_summary.csv
	echo "quality-filtered-reads," \$filt "("\$filt_perc"%)" >> ${sampleID}_summary.csv 
	"""
	
}

process umi_extraction {
        
	publishDir "$outDir/Data", mode: 'copy' //, pattern: '*_clusters_consensus.fasta'

        input:
        set sampleID, file(reads) from preprocessedReads
        val umiType from params.umiType

        output:
        tuple sampleID, file("*_extractedUMI.fasta") into extractedUmis

        script:
        if ( params.umiClustering )
        """
        extract_umis.py --max-error 3 --adapter-length 250 \
                --fwd-context "" \
                -o ${sampleID}_extractedUMI.fasta --tsv ${sampleID}_extractedUMI.tsv \
                $reads
        """
        else
        """
        touch dummy_extractedUMI.fasta
        """
}

process umi_clustering {
	
	publishDir "$outDir/Data", mode: 'copy', pattern: '*_clusters_consensus.fasta'

        input:
        set sampleID, file(umis) from extractedUmis
        val minlen from params.minlen
        val maxlen from params.maxlen
        val id from params.id

        output:
        tuple sampleID, file("${sampleID}_vsearch_clusters*") into clusters
        file("*_clusters_consensus.fasta") into consensus

        script:
        if ( params.umiClustering )
        """
        vsearch --clusterout_id --clusters ${sampleID}_vsearch_clusters \
                --centroids ${sampleID}_clusters_centroid.fasta \
                --consout ${sampleID}_clusters_consensus.fasta \
                --minseqlength $minlen --maxseqlength $maxlen --qmask none --threads 4 \
                --cluster_fast $umis --clusterout_sort --gapopen 0E/5I --gapext 0E/2I \
                --mismatch -8 --match 6 --iddef 0 --minwordmatches 0 --qmask none -id $id
        """
        else
        """
        touch ${sampleID}_vsearch_clusters_dummy
        touch dummy_clusters_consensus.fasta
        """
}

process umi_get_ubs {
        
	input:
        set sampleID, file(cluster) from clusters.transpose()

        output:
        tuple sampleID, env(size), file(cluster) into clusterTuples
        tuple sampleID, env(size), val(cluster.baseName) into getUBS

        script:
        if ( params.umiClustering )
        """
        size=`grep -c "^>" $cluster`
        """
        else
        """
        size=1000000
        """
}

ubsFile = outDir + "/ubs.txt"
cmd = "touch ${ubsFile}"
cmd.execute()
ubsFile = file(ubsFile)
ubsFile.append( getUBS.groupTuple( remainder:true ).collect().getVal() )


process umi_top_read {

        input:
        set sampleID, val(size), file(cluster) from clusterTuples

        output:
        tuple sampleID, file("*consensus.fasta") optional true into consensusUmis
        tuple sampleID, file("top.fasta"), file("*seqs.fasta") optional true into topUmis

        when:
        size.toInteger() >= params.ubs

        script:
        if ( size.toInteger() == 1 )
        """
        umi_to_seq.py -i $cluster -o ${cluster.baseName}_consensus.fasta
        """
        else if ( params.umiClustering )
        """
        vsearch --sortbysize ${cluster} --topn 1 --output top_umi.fasta
        umi_to_seq.py -i top_umi.fasta -o top.fasta
        umi_to_seq.py -i $cluster -o ${cluster.baseName}_seqs.fasta
        sed -i 's/>/>centroid_/' top.fasta
        """
        else
        """
        touch dummy_consensus.fasta
        """
}

process umi_polishing {

        container 'msanvicente/crisprgareq'
        errorStrategy 'ignore'

        input:
        set sampleID, file(top), file(seqs) from topUmis

        output:
        tuple sampleID, file(seqs), file("racon_2.fa") into raconPolished

        script:
        """
        minimap2 -x map-ont $top $seqs > ovlp_1.paf
        racon -t 4 -m 8 -x -6 -g -8 -w 500 --no-trimming $seqs ovlp_1.paf $top > racon_1.fa
        minimap2 -x map-ont racon_1.fa $seqs > ovlp_2.paf
        racon -t 4 -m 8 -x -6 -g -8 -w 500 --no-trimming $seqs ovlp_2.paf racon_1.fa > racon_2.fa
        """
}

process umi_consensus {

        errorStrategy 'ignore'

        input:
        set sampleID, file(seqs), file(racon) from raconPolished

        output:
        tuple sampleID, file("*consensus.fasta") into medakaConsensus

        script:
        """
        mini_align -i $seqs -r $racon -m -p map_1 -t 1
        medaka consensus map_1.bam consensus_1.hdf --threads 1 --model r941_min_high_g360 --chunk_len 6000
        mini_align -i $seqs -r $racon -m -p map_2 -t 1
        medaka consensus map_2.bam consensus_2.hdf --threads 1 --model r941_min_high_g360 --chunk_len 6000
        medaka stitch consensus_*.hdf ${seqs.baseName}_consensus.fasta
        """
}

process umi_final {
        
	publishDir "$outDir/Data", mode: 'copy', pattern: '*'
        beforeScript 'ulimit -Ss unlimited'

        input:
        set sampleID, file("*") from medakaConsensus.concat( consensusUmis ).groupTuple( remainder:true )

        output:
        tuple sampleID, file("${sampleID}_consensus.fasta") into finalConsensus

        script:
        """
        find . -type l -name '*consensus.fasta' -exec cat {} + > ${sampleID}_consensus.fasta
        """
}

process fasta_to_fastq {
	
	publishDir "$outDir/Data", mode: 'copy'
	
	input:
	set sampleID, file(consensusFasta) from finalConsensus

	output:
	tuple sampleID, file("*_consensus.fastq") into consensusFastq

	script:
	if ( params.umiClustering )
	"""
	fa2fq.py $consensusFasta > ${sampleID}_consensus.fastq
	"""
	else 
	"""
	touch dummy_consensus.fastq
	"""
}

simGEfiles = file( params.simGEfilesPath )

process simGE {
	
	publishDir "$outDir/Data", mode: 'copy'

	input:
	val seq from params.simGEsequence
	val cut from params.cutSite
	val sampleID from params.simGEid
	file sourceData from simGEfiles  

	output: 
	tuple sampleID, file("*_simulated.fastq") into simulated

	script:
	if( params.simGE )
	"""
	cp $simGEfiles/* .
	Rscript ./main.R $seq $cut 10000 ${sampleID}_simulated.fasta 
	fa2fq.py ${sampleID}_simulated.fasta > ${sampleID}_simulated.fastq	
	"""
	else
	"""
	touch dummy_simulated.fastq
	"""
}


process joinRead {
	
	publishDir "$outDir/Data", mode: 'copy'

	input:
	set sampleID, file(consensusReads), file(reads), sampleIDsim, file(simReads) from consensusFastq.join( noUMIpreprocessedReads ).combine( simulated )

	output:
	tuple sampleID, file("*_final-reads.fastq") into allProcessedReads
	tuple sampleID, file("*_final-reads.fastq") into readsSummary

	script:
	if ( params.umiClustering )
        """
	cp $consensusReads ${sampleID}_final-reads.fastq
	"""
	else if ( params.simGE )
	"""
	cp $simReads ${sampleID}_final-reads.fastq
	"""
	else
	"""
	cp $reads ${sampleID}_final-reads.fastq
	"""
}


process summary_reads{
	
	input:
	set sampleID, file(summary), file(finalReads) from summaryReport.join(readsSummary)
	
	output:
	set val(sampleID), file(summary) into summaryReportReads
	
	script:
	"""
	final_count=`expr \$(cat $finalReads | wc -l) / 4`
	echo "clustered-reads," \$final_count >> ${sampleID}_summary.csv
	"""
	
}


//index_human = file( params.indexHuman )
//index_mouse = file( params.indexMouse )

index_human = Channel.fromPath( params.indexHuman ).collect()
index_mouse = Channel.fromPath( params.indexMouse ).collect()


process selfReference {

        publishDir "$outDir/Data", mode: 'copy'

        input:
        set val(sampleID), file(reads), val(org), val(selfRefe), file(index_huma), file(index_mouse) from allProcessedReads.join( Channel.from(params.refOrganism) ).join( Channel.from(params.selfRef) ).combine( index_human ).combine( index_mouse )

        output:
        tuple sampleID, file(reads), file("*_reference.fasta") into refSequence

        script:
        if ( selfRefe )
        """        
        # Align reads against reference genome
        if [[ $org == human ]]; then
                bwa mem Human/hg38.fa.gz $reads > ${sampleID}.sam
        elif [[ $org == mouse ]]; then
                bwa mem Mouse/GCA_001632575.1_C3H_HeJ_v1.fa $reads > ${sampleID}.sam
        fi

	samtools view -Sb ${sampleID}.sam > ${sampleID}.bam
        samtools sort ${sampleID}.bam -o ${sampleID}.sorted.bam

        # Get coverage
        bedtools genomecov -bg -ibam ${sampleID}.sorted.bam > ${sampleID}.bed

        # Trim reads by minimal coverage (we skip it to avoid lossing samples with super low coverage)
        # trim_bed.py ${sampleID}.bed ${sampleID}_trimmed.bed #200

        # Merge adjacent sequences that have passed the minimal coverage
        # bedtools merge -i ${sampleID}_trimmed.bed > ${sampleID}_trimmed_merge.bed
        bedtools merge -i ${sampleID}.bed -c 4 -o mean > ${sampleID}_trimmed_merge.bed

        # Trim by longer unique interval
        # awk '{ printf ("%s\\t%s\\t%s\\t%s\\n", \$1,\$2,\$3,\$3-\$2) }' ${sampleID}_trimmed_merge.bed | sort -n -r -k 4 | head -1 > ${sampleID}_unique-longer.bed;

	# Trim by higher coverage
	sort -n -r -k 4 ${sampleID}_trimmed_merge.bed | head -1 > ${sampleID}_unique-longer.bed;

        # Get fasta sequence
        if [[ $org == human ]]; then
                bedtools getfasta -fi Human/hg38.fa -bed ${sampleID}_unique-longer.bed -fo ${sampleID}_reference.fasta
        elif [[ $org == mouse ]]; then
                bedtools getfasta -fi Mouse/GCA_001632575.1_C3H_HeJ_v1.fa -bed ${sampleID}_unique-longer.bed -fo ${sampleID}_reference.fasta
        fi

	#change the title of reference
        sed -i "1s/.*/>${sampleID}/" ${sampleID}_reference.fasta
        
        """
        else
        """
        touch dummy_reference.fasta
        """

}


process givenRefSeq {

        publishDir "$outDir/Data", mode: 'copy'

        input:
        set sampleID, file(refSeq), val(selfRefe) from Channel.fromPath(params.referenceFasta).map{ file -> tuple(file.baseName, file) }.join( Channel.from(params.selfRef) )

        output:
        tuple sampleID, file("*_reference.fasta") into givenRefSequence

        script:
        if ( ! selfRefe )
        """
        cp $refSeq ${sampleID}_reference.fasta
	#change the title of reference
        sed -i "1s/.*/>${sampleID}/" ${sampleID}_reference.fasta
        """
        else
        """
        touch dummy_reference.fasta
        """

}

process givenTemplate {

        publishDir "$outDir/Data", mode: 'copy'

        input:
        set sampleID, file(templateSeq) from Channel.fromPath(params.template_seq).map{ file -> tuple(file.baseName, file) }

        output:
        tuple sampleID, file("*_template.fasta") into templateSequence

        script:
        if ( params.template )
        """
        cp $templateSeq ${sampleID}_template.fasta
        """
        else
        """
        touch dummy_template.fasta
        """

}

process gRNAoriantationRef {
	
	publishDir "$outDir/Data", mode: 'copy'

	input:
        set val(sampleID), file(reads), file(refSeq), file(givRefSeq), val(gRNAseq) from refSequence.join( givenRefSequence ).join( Channel.from(params.gRNAseq) )

	output:
	tuple sampleID, file(reads), file("*_reference-correctOrient.fasta") into refSeqOriented

	script:
	"""
	# Check if gRNA is in reference sequence (same orientation between gRNA and reference) and return reverse complement if necessary
	revComp-fasta.R ${sampleID}_reference.fasta ${sampleID}_reference-correctOrient.fasta $gRNAseq;
	"""

}

process alignment {
        
	publishDir "$outDir/Data", mode: 'copy'

        input:
        set sampleID, file(reads), file(refSeq) from refSeqOriented
	val alignMethod from params.aligner

        output:
        tuple sampleID, file(refSeq), file(reads), file("*.sorted.bam"), file("*.sorted.bam.bai") into align
        tuple sampleID, file("*.sorted.bam") into summaryAlignment
        file("*.fasta.fai") into genomeBrowser

	script:
	if ( alignMethod == "minimap2" )
		"""
		minimap2 -d ${sampleID}_reference.mmi ${refSeq};
		minimap2 -a -A 29 -B 17 -O 25 -E 2 ${sampleID}_reference.mmi ${reads} > ${sampleID}.sam;
		samtools view -bS ${sampleID}.sam -o ${sampleID}.bam;
		samtools sort ${sampleID}.bam -o ${sampleID}.sorted.bam;
		samtools index ${sampleID}.sorted.bam;
        	samtools faidx ${refSeq} -o ${refSeq}.fai;
		"""

	else if ( alignMethod == "bwa" )
		"""
		bwa index ${refSeq};
        	bwa mem ${refSeq} ${reads} > ${sampleID}.sam;
		samtools view -bS ${sampleID}.sam -o ${sampleID}.bam;
		samtools sort ${sampleID}.bam -o ${sampleID}.sorted.bam;
		samtools index ${sampleID}.sorted.bam;
        	samtools faidx ${refSeq} -o ${refSeq}.fai;
		"""

	else if ( alignMethod == "bowtie2" )
		"""
		bowtie2-build -f ${refSeq} ${sampleID}_reference 
		bowtie2 -x ${sampleID}_reference -U ${reads} -S ${sampleID}.sam;
		samtools view -bS ${sampleID}.sam -o ${sampleID}.bam; 
		samtools sort ${sampleID}.bam -o ${sampleID}.sorted.bam;
		samtools index ${sampleID}.sorted.bam;
        	samtools faidx ${refSeq} -o ${refSeq}.fai;
		"""
}

process summary_aligned{
	
	publishDir "$outDir/Data", mode: 'copy'
	
	input:
	set sampleID, file(alignedReads), file(summary) from summaryAlignment.join(summaryReportReads)
	
	output:
	tuple sampleID, file(summary) into finalSummary
	
	script:
	"""
	mapped_count=`samtools view -c -b -F 4 $alignedReads`
	total=`samtools view -c -b $alignedReads` 
	perc=`echo "scale=1;(\$mapped_count/\$total)*100" | bc`
	echo "aligned-reads," \$mapped_count "("\$perc"%)" >> ${sampleID}_summary.csv 
	"""
	
}


process variantCalling {
	
	publishDir "$outDir/Results", mode: 'copy'
	
	input:
	set sampleID, file(refSeq), file(reads), file(alignedReads), file(indexFile), val(gRNA), file(summary), file(template), val(relCutSite) from align.join(Channel.from(params.gRNAseq)).join(finalSummary).join(templateSequence).join(Channel.from(params.protCutSite))
        val spike from params.spike
        val mock from params.mock
	
	output:
	tuple sampleID, file(refSeq), val(gRNA), file("${sampleID}_*indels.csv"), file("${sampleID}_subs-perc.csv"), val(relCutSite) into indels
	file("${sampleID}_*indels.csv") into indels2correct
	tuple file("${sampleID}*.html"), file("*edits.csv") into edition //file("${sampleID}_edit-perc.csv"), 
	    file("*cutSite.json") into cutForGenomeBrowser

	script:
	"""
	VC_parser-cigar.R ${alignedReads} ${sampleID} ${refSeq} ${gRNA} ${sampleID} ${template} ${spike} ${summary} ${relCutSite} ${mock};
	"""

}

test = indels2correct.map{denoise -> [denoise.simpleName, denoise]}.groupTuple()


process noiseCleaning {
	
	publishDir "$outDir/Results", mode: 'copy'
	
	input:
	//set sampleID, file(boo) from test
	val mock from params.mock
	set sampleID, file(boo) from test
	
	output:
	file("*CorrectedIndels.csv") into Newindels	
	file("*MockSummary.csv") into MockSummary
	file("*Mockedits.csv") into Mockthings2
	file("*MvsTdel.html") into Mockthings3
	file("*MvsTins.html") into Mockthings
	
	when: 
	mock==true
	
	script:
	"""
	noisemock.R $boo;
	"""
}

process getSubstractionPlots {
	
	publishDir "$outDir/Results", mode: 'copy'

	input:
	val mock from params.mock
	file(corrected) from Newindels

	output:
	file("*_insSubstraction.html") into SubstractionPlot
    file("*_delsSubstraction.html") into SubstractionPlot2
	
	when: 
	mock==true
	
	script:
	"""
	SubstractionPlot.R $corrected
	"""
}

/*
process getSubstitutions {

	publishDir "$outDir/Results", mode: 'copy'

	input:
	set sampleID, file(refSeq), val(gRNA), file(reads), file(indels), file(subs_perc) from indels

	output: 
	tuple sampleID, file(refSeq), val(gRNA), file(indels),  file(subs_perc), file("substitutions_info.RData") into allEdits

	script:
	"""
	fq2fa.R $reads read_1
	blat -t=dna -q=dna -minIdentity=90 -minScore=10 -out=pslx read_1.fasta $refSeq read_1.psl
	substitution-VC.R $refSeq $gRNA	
	"""
}
*/


process getSamplePlots {
	
	publishDir "$outDir/Results", mode: 'copy'

	input:
	set sampleID, file(refseq), val(gRNA), file(editsInfo), file(substitutions_pileup), val(relCutSite) from indels //allEdits

	output:
	tuple sampleID, file("*_Deletions.html"), file("*_Insertions.html"), file("*_top-alleles_LOGO.png"), file("*_top.html"), file("*_subs-perc_plot.png"), file("*_subs-perc_plot_LOGO.png") into samplePlots
	file(editsInfo) into editsCSV
	file("*.json") into dynamicTable
	script:
	"""
	get_plots.R $sampleID $editsInfo $refseq $gRNA $substitutions_pileup $relCutSite
	"""
}



process samplesPreHeatMap {
	
	publishDir "$outDir/Results", mode: 'copy'

	input:
	file(editsInfo) from editsCSV.collect( )

	output:
	tuple file("all-samples_indels.csv"), env(count) into preHeatMap
	
	script:
	"""
	count=0
	for sample in ${editsInfo}
	do
		let count+=1
		cat \$sample >> all-samples_indels.csv
	done
	"""
}


sampleNames = file( params.samplesNames_file )

process samplesDistanceHeatmap {
	
	publishDir "$outDir/Results", mode: 'copy'

	input:
	tuple file(samplesIndels), val(count) from preHeatMap
	file samples_names from sampleNames

	output:
	file("*.html") into heatMaps

	script:
	if ( count != "1" && count != "2" )
	"""
	heatmap_samplesDist.R $samplesIndels $samples_names
	"""
	else
	"""
	touch dummy.html
	"""
}
