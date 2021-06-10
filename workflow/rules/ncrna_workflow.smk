from collections import defaultdict

include: "bowtie2.smk"
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4


#directories, files, = glob_wildcards("analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam")

directories, files, = glob_wildcards("data/raw/{samplegroup}/{sample}.fastq.gz")

gene=["gencode","piRNA","miRNA_precursor","miRNA","tRNA"]


configfile: "config/config.yaml"


file_dict=defaultdict(list)

directory_workdir=os.getcwd()


#assign files to correct directories
for i in zip(directories,files):
	file_dict[i[0]].append(i[1])

	

rule all_ncrna:
	input:
		expand("results/count_tables/{samplegroup}/{gene}.tsv",samplegroup=set(directories),gene=gene),
		expand(["results/count_tables/{samplegroup}/" + x + ".tsv" for x in ["lincRNA","misc_RNA","protein_coding","snoRNA","snRNA","scaRNA","antisense"]],samplegroup=set(directories))




rule exon_table_from_gencode:
	input: 
			bam="analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam",
			gencode=config["gencode"]

	output:
		"analyses/bowtie_mappings_genome_multi/txt_tables_per_file/{samplegroup}/" + x + "/{sample}." + x + ".txt" for x in ["gencode"]

	conda:
		"envs/main.yaml"
	shell:
		'''	
		input_name=$(basename {input.bam})
		samplename=${{input_name/.sorted.bam/}}
		workflow/scripts/rnaseq_tool.py {input.bam} -g {input.gencode} --collapsed --filter -f exon -i gene_name | awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' > {output}
		'''

rule pirna_from_pirbase:
	input:
		bam="analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam",
		pirbase=config["pirbase"]
	output:
		"analyses/bowtie_mappings_genome_multi/txt_tables_per_file/{samplegroup}/" + x + "/{sample}." + x + ".txt" for x in ["piRNA"]
	shell:
		'''
        input_name=$(basename {input.bam})
		samplename=${{input_name/.sorted.bam/}}
			workflow/scripts/rnaseq_tool.py {input.bam} -g {input.pirbase} --collapsed --filter -f piRNA -i name | awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' > {output}

		'''

rule mirna_and_precursor_from_mirBase:
	input:
		bam="analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam",
		mirbase=config["mirbase"]
	output:
		"analyses/bowtie_mappings_genome_multi/txt_tables_per_file/{samplegroup}/" + x + "/{sample}." + x + ".txt" for x in ["miRNA","miRNA_precursor"]
	shell:  
                '''
                	input_name=$(basename {input.bam})
		samplename=${{input_name/.sorted.bam/}}

		workflow/scripts/rnaseq_tool.py {input.bam} -g {input.mirbase} --collapsed --filter  -f miRNA -i Name |  awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' > {output[0]}
		workflow/scripts/rnaseq_tool.py {input.bam} -g {input.mirbase} --collapsed --filter  -f miRNA_primary_transcript -i Name |  awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' > {output[1]}

                '''




rule trna_from_gencode:
	input:
		bam="analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam",
		trnagencode=config["trnagencode"]
	output:
		"analyses/bowtie_mappings_genome_multi/txt_tables_per_file/{samplegroup}/" + x + "/{sample}." + x + ".txt" for x in ["tRNA"]
	
	shell:
		'''
                	input_name=$(basename {input.bam})
		samplename=${{input_name/.sorted.bam/}}

		workflow/scripts/rnaseq_tool.py {input.bam} -g {input.trnagencode} --collapsed --filter -f tRNA -i gene_name | awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' > {output};

                '''

rule create_final_count_tables:
	input:
		lambda wildcards : ["analyses/bowtie_mappings_genome_multi/txt_tables_per_file/" + wildcards.samplegroup + "/" + wildcards.gene + "/" + x + "." + wildcards.gene + ".txt" for x in file_dict[wildcards.samplegroup]]

	output:
		"results/count_tables/{samplegroup}/{gene}.tsv"

	shell:
		'''
		workflow/scripts/generic_table_creator.R {wildcards.gene} analyses/bowtie_mappings_genome_multi/txt_tables_per_file/{wildcards.samplegroup}/{wildcards.gene}
		mv analyses/bowtie_mappings_genome_multi/txt_tables_per_file/{wildcards.samplegroup}/{wildcards.gene}/{wildcards.gene}.tsv results/count_tables/{wildcards.samplegroup}/{wildcards.gene}.tsv
		'''


rule separate_gene_types:
	input:
		"results/count_tables/{samplegroup}/gencode.tsv"
	output:
		["results/count_tables/{samplegroup}/" + x + ".tsv" for x in ["lincRNA","misc_RNA","protein_coding","snoRNA","snRNA","scaRNA","antisense"]]
	run:
		shell("""
		cd results/count_tables/{wildcards.samplegroup}/
		{directory_workdir}/workflow/scripts/separate_gene_types.R gencode.tsv

		""")

