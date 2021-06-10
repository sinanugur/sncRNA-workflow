directories, files, = glob_wildcards("data/raw/{samplegroup}/{sample}.fastq.gz")

humangenome="/home/sium/data/humangenome/hg38"

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

rule all_bowtie2:
	input:
		expand("analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam",zip,samplegroup=directories,sample=files),
		"results/file_statistics.csv",
		expand("data/stats/{samplegroup}/{sample}.trimmed.FastQC",zip,samplegroup=directories,sample=files)

rule adapter_trimming_fastp:
	input:
		"data/raw/{samplegroup}/{sample}.fastq.gz"

	output:
		"data/trimmed/{samplegroup}/{sample}.trimmed.fastq.gz"
	
	threads: 15

	shell:
		"""

		fastp -i {input} -o {output} -w {threads}

		"""

rule collapsing_reads:
	input:
		"data/trimmed/{samplegroup}/{sample}.trimmed.fastq.gz"
	
	output:
		"data/collapsed/{samplegroup}/{sample}.trimmed.collapsed.fasta.gz"

	shell:
		"""
		
		zcat {input} | fastx_collapser | gzip > {output}
				
		"""

rule bowtie2_mapping:
	input:
		"data/collapsed/{samplegroup}/{sample}.trimmed.collapsed.fasta.gz"

	output:
		"analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam"

	threads: 15

	shell:
		#"bowtie_map_local_all.sh /home/sium/data/humangenome/hg38 {input} {output}"
		"""
		bowtie2 --sensitive-local -k 10 -f -p {threads} -x {humangenome} -U <(zcat {input}) | samtools view -bS - | samtools sort - -o {output}
		"""

rule file_stats:
	input:
		"data/raw/{samplegroup}/{sample}.fastq.gz",
		"data/trimmed/{samplegroup}/{sample}.trimmed.fastq.gz",
		"data/collapsed/{samplegroup}/{sample}.trimmed.collapsed.fasta.gz"

	output:
		"data/stats/{samplegroup}/{sample}.sample.stats.txt"

	shell:
		"""
		md5=$(md5sum {input[0]} | awk '{{print $1}}')
		seq_count=$(zcat {input[0]} | echo $((`wc -l`/4)))
		trimmed_seq_count=$(zcat {input[1]} | echo $((`wc -l`/4)))
		collapsed_reads=$(zcat {input[2]} | grep ">" -c)
		
	
		filename=$(basename {input[0]})
		gc=$(seqtk fqchk {input[0]} | awk '/ALL/{{print $4+$5}}')
				
		filesize=$(du -h --block-size=M $physical | cut -f1)

		echo -e "md5sum\tFilename\t#Reads\t#Reads After Trimming\t#Collapsed Reads\t#File Size(MB)\tGC" > {output} 
		echo -e $md5"\t"$filename"\t"$seq_count"\t"$trimmed_seq_count"\t"$collapsed_reads"\t"$filesize"\t"$gc >> {output}

		"""

rule fastqc_stats:
	input:
		"data/trimmed/{samplegroup}/{sample}.trimmed.fastq.gz"

	output:
		directory("data/stats/{samplegroup}/{sample}.trimmed.FastQC")

	threads: 15

	shell:
		"""
		mkdir {output};
		fastqc --nogroup --extract -o {output} -t {threads} {input};
		"""

rule combine_file_stats:
	input:
		expand("data/stats/{samplegroup}/{sample}.sample.stats.txt",zip,samplegroup=directories,sample=files)

	output:
		"results/file_statistics.csv"
	
	run:
		shell("""echo "md5sum\tphysical\tfilename\tcase_control\tsample\ttotal_reads\ttotal_reads_after_trimming\tcollapsed_reads\tspike_in_counts\tfile_size(MB)\tgc" > {output};""")	
		for i in input:
			shell("""cat {i} | grep -v "#" >> {output};""")
