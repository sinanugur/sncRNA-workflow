from collections import defaultdict
import os.path


if config.get('directory') is None:
	directories, files, = glob_wildcards("data/raw/{samplegroup}/{sample}.fastq.gz")
	#glob_wildcards("WGS/tests/fastq/{group}/{sample}.fastq.gzadr_R1.fastq")
else: #a directory name can be provided for a single group run
	files, = glob_wildcards("data/raw/" + config.get('directory') + "/{sample}.fastq.gz")
	directories=[config.get('directory')]*len(files)



mirbase_directory="/WORKING/databases/miRBase/"


print(directories)
print(files)

rule all:
	input:
		expand("analyses/mirtop/{samplegroup}/bowtie/{sample}.mature.bam",zip,samplegroup=directories,sample=files)
	

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


rule mirtop_collapse_or_unzip:
	input:
		"data/trimmed/{samplegroup}/{sample}.trimmed.fastq.gz"
	
	output:
		"analyses/mirtop/{samplegroup}/collapse/{sample}.trimmed.fastq"

	params:
		o="analyses/mirtop/{samplegroup}/collapse"
	
	shell:
		"""
	
		seqcluster collapse -f {input}  -m 1 --min_size 17 -o {params.o}

		"""


rule bowtie_mature:
	input:
		"data/trimmed/{samplegroup}/{sample}.trimmed.fastq.gz"
		
	output:
		unmapped="analyses/mirtop/{samplegroup}/bowtie/{sample}.mature_unmapped.fq",
		bam="analyses/mirtop/{samplegroup}/bowtie/{sample}.mature.bam"

	params:
		mature_base= "/WORKING/databases/miRBase/hsa-mature"
	threads: 15
	shell:
		"""
		    bowtie \\
        {params.mature_base} \\
        -q <(zcat {input}) \\
        -p {threads} \\
        -t \\
        -k 50 \\
        --best \\
        --strata \\
        -e 99999 \\
        --chunkmbs 2048 \\
        --un {output.unmapped}\\
        -S  \\
        | samtools view -bS - > {output.bam}


		"""


rule bowtie_hairpin:
	input:
		"analyses/mirtop/{samplegroup}/bowtie/{sample}.mature_unmapped.fq"
		
	output:
		unmapped="analyses/mirtop/{samplegroup}/bowtie/{sample}.hairpin_unmapped.fq",
		bam="analyses/mirtop/{samplegroup}/bowtie/{sample}.hairpin.bam"

	params:
		hairpin_base= "/WORKING/databases/miRBase/hsa-hairpin"
	threads: 15
	shell:
		"""
		    bowtie \\
        {params.hairpin_base} \\
        -q <(zcat {input}) \\
        -p {threads} \\
        -t \\
        -a \\
        --best \\
        --strata \\
        -e 99999 \\
        --chunkmbs 2048 \\
        --un {output.unmapped}\\
        -S  \\
        | samtools view -bS - > {output.bam}


		"""


rule bowtie_hairpin_collapsed:
	input:
		"analyses/mirtop/{samplegroup}/collapse/{sample}.trimmed.fastq"
		
	output:
		bam="analyses/mirtop/{samplegroup}/bowtie/{sample}.hairpin_collapsed.bam"

	params:
		hairpin_base= "/WORKING/databases/miRBase/hsa-hairpin"
	threads: 15
	shell:
		"""

    bowtie \\
        {params.hairpin_base} \\
        -p {threads} \\
        -t \\
        -k 50 \\
        -a \\
        --best \\
        --strata \\
        -e 99999 \\
        --chunkmbs 2048 \\
        -q <(cat {input}) \\
        -S \\
        | samtools view -bS - > {output.bam}
    """

		"""