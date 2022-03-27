from collections import defaultdict

include: "bowtie2.smk"
include: "bowtie1.smk"
include: "prepare_databases.smk"
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4


#directories, files, = glob_wildcards("analyses/bowtie_mappings_genome_multi/{sample}.sorted.bam")

#files, = glob_wildcards("data/{sample}.fastq.gz")

#gene=["gencode","piRNA","miRNA_precursor","miRNA","tRNA"]


configfile: "config/config.yaml"


file_dict=defaultdict(list)

directory_workdir=os.getcwd()


#assign files to correct directories
#for i in zip(directories,files):
#	file_dict[i[0]].append(i[1])

	

rule all_ncrna:
	input:
		#expand("results/count_tables/{gene}.tsv",gene=gene),
		expand(["results/count_tables/" + x + ".tsv" for x in ["lincRNA","misc_RNA","protein_coding","snoRNA","snRNA","scaRNA","antisense"]])




rule exon_table_from_gencode:
	input: 
		bam="analyses/bowtie_mappings_genome_multi/{sample}.sorted.bam",
		gencode="databases/gencode.gff3"

	output:
		"analyses/bowtie_mappings_genome_multi/txt_tables_per_file/" + x + "/{sample}." + x + ".txt" for x in ["gencode"]

	conda:
		"../envs/main.yaml"
	shell:
		'''	
		input_name=$(basename {input.bam})
		samplename=${{input_name/.sorted.bam/}}
		workflow/scripts/rnaseq_tool.py {input.bam} -g {input.gencode} --collapsed --filter -f exon -i gene_name | awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' > {output}
		'''

rule pirna_from_pirbase:
	input:
		bam="analyses/bowtie_mappings_genome_multi/{sample}.sorted.bam",
		pirbase="databases/piRNA.gtf"
	output:
		"analyses/bowtie_mappings_genome_multi/txt_tables_per_file/" + x + "/{sample}." + x + ".txt" for x in ["piRNA"]

	conda:
		"../envs/main.yaml"
	shell:
		'''
        input_name=$(basename {input.bam})
		samplename=${{input_name/.sorted.bam/}}
			workflow/scripts/rnaseq_tool.py {input.bam} -g {input.pirbase} --collapsed --filter -f piRNA -i name | awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' > {output}

		'''

rule mirna_and_precursor_from_mirBase:
	input:
		bam="analyses/bowtie_mappings_genome_multi/{sample}.sorted.bam",
		mirbase=config["mirbase"]
	output:
		"analyses/bowtie_mappings_genome_multi/txt_tables_per_file/" + x + "/{sample}." + x + ".txt" for x in ["miRNA","miRNA_precursor"]


	conda:
		"../envs/main.yaml"
	shell:  
                '''
                	input_name=$(basename {input.bam})
		samplename=${{input_name/.sorted.bam/}}

		workflow/scripts/rnaseq_tool.py {input.bam} -g {input.mirbase} --collapsed --filter  -f miRNA -i Name |  awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' > {output[0]}
		workflow/scripts/rnaseq_tool.py {input.bam} -g {input.mirbase} --collapsed --filter  -f miRNA_primary_transcript -i Name |  awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' > {output[1]}

                '''




rule trna_from_gencode:
	input:
		bam="analyses/bowtie_mappings_genome_multi/{sample}.sorted.bam",
		trnagencode=config["trnagencode"]
	output:
		"analyses/bowtie_mappings_genome_multi/txt_tables_per_file/" + x + "/{sample}." + x + ".txt" for x in ["tRNA"]

	conda:
		"../envs/main.yaml"
	
	shell:
		'''
                	input_name=$(basename {input.bam})
		samplename=${{input_name/.sorted.bam/}}

		workflow/scripts/rnaseq_tool.py {input.bam} -g {input.trnagencode} --collapsed --filter -f tRNA -i gene_id | awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' > {output};

                '''

rule create_final_count_tables:
	input:
		lambda wildcards : ["analyses/bowtie_mappings_genome_multi/txt_tables_per_file/" + wildcards.gene + "/" + s + "." + wildcards.gene + ".txt" for s in files]

	output:
		"results/count_tables/{gene}.tsv"

	conda:
		"../envs/main.yaml"

	shell:
		'''
		#Rscript --no-environ workflow/scripts/generic_table_creator.R {wildcards.gene} analyses/bowtie_mappings_genome_multi/txt_tables_per_file/{wildcards.gene}
		#mv analyses/bowtie_mappings_genome_multi/txt_tables_per_file/{wildcards.gene}/{wildcards.gene}.tsv results/count_tables/{wildcards.gene}.tsv

		echo "{output} {input}" | xargs Rscript --no-environ ./workflow/scripts/sncrna_table_creator.R
		'''


rule separate_gene_types:
	input:
		"results/count_tables/gencode.tsv",
		"databases/gencode.gene-name.csv"
	output:
		["results/count_tables/" + x + ".tsv" for x in ["lncRNA","misc_RNA","protein_coding","snoRNA","snRNA","scaRNA"]]


	conda:
		"../envs/main.yaml"
	shell:
		"""
		cd results/count_tables/
		Rscript --no-environ {directory_workdir}/workflow/scripts/separate_gene_types.R gencode.tsv

		"""

