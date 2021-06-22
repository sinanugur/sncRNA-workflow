#directories, files, = glob_wildcards("data/raw/{sample}.fastq.gz")


# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

configfile: "config/config.yaml"

humangenome=config["humangenome"]

rule all_bowtie2:
	input:
		expand("analyses/bowtie_mappings_genome_multi/{sample}.sorted.bam",sample=files),
		"results/file_statistics.csv",
		expand("analyses/stats/{sample}.trimmed.FastQC",sample=files)

rule adapter_trimming_fastp:
	input:
		"data/{sample}.fastq.gz"

	output:
		"analyses/trimmed/{sample}.trimmed.fastq.gz"

	conda:
		"../envs/main.yaml"
	
	threads: 15

	shell:
		"""

		fastp -i {input} -o {output} -w {threads}

		"""

rule collapsing_reads:
	input:
		"analyses/trimmed/{sample}.trimmed.fastq.gz"
	
	output:
		"analyses/collapsed/{sample}.trimmed.collapsed.fasta.gz"

	conda:
		"../envs/main.yaml"

	shell:
		"""
		
		zcat {input} | fastx_collapser | gzip > {output}
				
		"""

rule bowtie2_mapping:
	input:
		"analyses/collapsed/{sample}.trimmed.collapsed.fasta.gz",
		multiext(
            "databases/hg38",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        )

	output:
		"analyses/bowtie_mappings_genome_multi/{sample}.sorted.bam"

	conda:
		"../envs/main.yaml"

	threads: 15

	shell:
		"""
		bowtie2 --sensitive-local -k 10 -f -p {threads} -x {humangenome} -U <(zcat {input[0]}) | samtools view -bS - | samtools sort - -o {output}
		"""

rule file_stats:
	input:
		"data/{sample}.fastq.gz",
		"analyses/trimmed/{sample}.trimmed.fastq.gz",
		"analyses/collapsed/{sample}.trimmed.collapsed.fasta.gz"

	output:
		"analyses/stats/{sample}.sample.stats.txt"

	conda:
		"../envs/main.yaml"

	shell:
		"""
		md5=$(md5sum {input[0]} | awk '{{print $1}}')
		seq_count=$(zcat {input[0]} | echo $((`wc -l`/4)))
		trimmed_seq_count=$(zcat {input[1]} | echo $((`wc -l`/4)))
		collapsed_reads=$(zcat {input[2]} | grep ">" -c)
		if [ -L {input[0]} ]; then physical=$(readlink {input[0]}); else physical={input[0]}; fi
		
	
		filename=$(basename {input[0]})
		gc=$(seqtk fqchk {input[0]} | awk '/ALL/{{print $4+$5}}')
				
		filesize=$(du -h --block-size=M $physical | cut -f1)

		echo -e "md5sum\tFilename\t#Reads\t#Reads After Trimming\t#Collapsed Reads\t#File Size(MB)\tGC" > {output} 
		echo -e $md5"\t"$filename"\t"$seq_count"\t"$trimmed_seq_count"\t"$collapsed_reads"\t"$filesize"\t"$gc >> {output}

		"""

rule fastqc_stats:
	input:
		"analyses/trimmed/{sample}.trimmed.fastq.gz"

	output:
		directory("analyses/stats/{sample}.trimmed.FastQC")

	conda:
		"../envs/main.yaml"

	threads: 15

	shell:
		"""
		mkdir {output};
		fastqc --nogroup --extract -o {output} -t {threads} {input};
		"""

rule combine_file_stats:
	input:
		expand("analyses/stats/{sample}.sample.stats.txt",sample=files)

	output:
		"results/file_statistics.csv"
	
	run:
		shell("""echo "md5sum\tphysical\tfilename\tcase_control\tsample\ttotal_reads\ttotal_reads_after_trimming\tcollapsed_reads\tspike_in_counts\tfile_size(MB)\tgc" > {output};""")	
		for i in input:
			shell("""cat {i} | grep -v "#" >> {output};""")
