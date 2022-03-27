'''
@author: Sinan U. Umu, PhD, sinanugur@gmail.com
'''

rule bowtie1_mapping:
	input:
		"analyses/trimmed/{sample}.trimmed.fastq.gz",
		multiext(
			"databases/hg38",
			".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt",
		)

	output:
		"analyses/bowtie1_mappings/{sample}.sorted.bam"
	threads: 15 
	shell:
		"""
		bowtie -k 10 --best --strata --chunkmbs 2000 -p {threads} --sam {humangenome} -q <( zcat {input[0]} ) | samtools view -bS - | samtools sort - -o {output}
		"""


