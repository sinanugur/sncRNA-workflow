from collections import defaultdict

include: "bowtie2.smk"
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4


#directories, files, = glob_wildcards("analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam")
directories, files, = glob_wildcards("data/raw/{samplegroup}/{sample}.fastq.gz")
gene=["gencode.uniq.exon","pirbase.uniq","tRNAgencode.uniq"]


from collections import defaultdict

include: "bowtie2.smk"
include: "prepare_databases.smk"
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4


files, = glob_wildcards("data/{sample}.fastq.gz")

configfile: "config/config.yaml"


rule all_statistics:
	input:
		expand("results/statistics/{samplegroup}/{gene}.count.table.csv",samplegroup=set(directories),gene=gene),




rule feature_counts:
	input:
		bams=expand("analyses/bowtie1_mappings/{sample}.sorted.bam",sample=files),
		gencode="databases/gencode.gff3"
	output:
		"analyses/featurecounts/gencode.txt"
	threads: 20
	shell:
		"""
		featureCounts -O -g gene_type -T {threads} -s 1 -a /WORKING/databases/gencode/gencode.v26.annotation.gff3 -o {output}  {input.bams}
		"""

