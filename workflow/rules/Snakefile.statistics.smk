from collections import defaultdict


# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4


directories, files, = glob_wildcards("analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam")
gene=["gencode.uniq.exon","pirbase.uniq","tRNAgencode.uniq"]


configfile: "workflows/config.yaml"
version=str(config["version"])

directory_bowtie2="analyses/bowtie_mappings_genome_multi/count_tables/" + version + "/"
directory_seqbuster="analyses/seqbuster/count_tables/"
directory_results="results/" + version + "/"
directory_trf="analyses/MINTmap/count_tables/"

file_dict=defaultdict(list)

#assign files to correct directories
for i in zip(directories,files):
	file_dict[i[0]].append(i[1])

	

rule all:
	input:
		expand(directory_bowtie2 + "{samplegroup}/{gene}.count.table.csv",samplegroup=set(directories),gene=gene),




rule uniq_mapping_gencode:
	input:
		bam="analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam",
		gencode=config["gencode"]
	output:
		"analyses/bowtie_mappings_genome_multi/txt_tables_per_file/" + version + "/{samplegroup}/" + g +  "/{sample}." + g + ".txt" for g in ["gencode.uniq.exon"]

	params:
		par="--collapsed --unique -f exon -i gene_type"
	shell:
		"""
		samplename=$(bin/sample_name.sh {input})
		bin/rnaseq_tool.py {input.bam} -g {input.gencode} {params.par} | awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' >  {output}
		"""

rule uniq_mapping_trna:
	input:
		bam="analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam",
		trnagencode=config["trnagencode"].strip("gff3") + "forstatistics.gff3"

	output:
		"analyses/bowtie_mappings_genome_multi/txt_tables_per_file/" + version +  "/{samplegroup}/" + g +  "/{sample}." + g + ".txt" for g in ["tRNAgencode.uniq"]	
	params:
		par = "--collapsed --unique -f tRNA -i stat_type"
	
	shell:
		"""
                samplename=$(bin/sample_name.sh {input})
                bin/rnaseq_tool.py {input.bam} -g {input.trnagencode} {params.par} | awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' >  {output}
		"""


rule uniq_mapping_pirna:
	input:
		bam="analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam",
		pirbase=config["pirbase"]
	output:
		"analyses/bowtie_mappings_genome_multi/txt_tables_per_file/" + version + "/{samplegroup}/" + g +  "/{sample}." + g + ".txt" for g in ["pirbase.uniq"]	
	params:
		par = "--collapsed --unique -f piRNA -i gene_type"

	shell:
		"""
		samplename=$(bin/sample_name.sh {input})
		bin/rnaseq_tool.py {input.bam} -g {input.pirbase} {params.par} | awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' >  {output}
        """


rule create_final_uniq_count_tables:
	input:
		lambda wildcards : ["analyses/bowtie_mappings_genome_multi/txt_tables_per_file/" + wildcards.version + "/" + wildcards.samplegroup + "/" + wildcards.gene + "/" + x + "." + wildcards.gene + ".txt" for x in file_dict[wildcards.samplegroup]]

	output:
		 "analyses/bowtie_mappings_genome_multi/count_tables/{version}/{samplegroup}/{gene}.count.table.csv"

	shell:
		'''
		bin/generic_table_creator.R {wildcards.gene} analyses/bowtie_mappings_genome_multi/txt_tables_per_file/{wildcards.version}/{wildcards.samplegroup}/{wildcards.gene}
		mv analyses/bowtie_mappings_genome_multi/txt_tables_per_file/{wildcards.version}/{wildcards.samplegroup}/{wildcards.gene}/{wildcards.gene}.count.table.csv analyses/bowtie_mappings_genome_multi/count_tables/{wildcards.version}/{wildcards.samplegroup}/{wildcards.gene}.count.table.csv
		'''
		

