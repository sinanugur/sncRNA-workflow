from collections import defaultdict


# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4


directories, files, = glob_wildcards("analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam")
gene=["exons.gencode","gene.gencode","pirbase","miRBaseprecursor","miRBasemiRNA","mirgenedbmirna","mirgenedbprecursor","tRNAgencode"]


configfile: "workflows/config.yaml"
version=str(config["version"])

file_dict=defaultdict(list)

#assign files to correct directories
for i in zip(directories,files):
	file_dict[i[0]].append(i[1])

	

rule all:
	input:
		expand("analyses/bowtie_mappings_genome_multi/count_tables/" + version + "/{samplegroup}/{gene}.count.table.csv",samplegroup=set(directories),gene=gene)

rule exon_table_from_gencode:
	input: 
			bam="analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam",
			gencode=config["gencode"]

	output:
		"analyses/bowtie_mappings_genome_multi/txt_tables_per_file/" + version + "/{samplegroup}/" + x + "/{sample}." + x + ".txt" for x in ["exons.gencode"]
	shell:
		'''	
		samplename=$(bin/sample_name.sh {input.bam})
		bin/rnaseq_tool.py {input.bam} -g {input.gencode} --collapsed --filter -f exon -i gene_name | awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' > {output}
		'''


rule gene_table_from_gencode:
	input:
		bam="analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam",
		gencode=config["gencode"]

	output:
		"analyses/bowtie_mappings_genome_multi/txt_tables_per_file/" + version + "/{samplegroup}/" + x + "/{sample}." + x + ".txt" for x in ["gene.gencode"]

	shell:
		'''
                samplename=$(bin/sample_name.sh {input.bam})
                bin/rnaseq_tool.py {input.bam} -g {input.gencode} --collapsed --filter -f gene -i gene_name | awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' > {output}
                '''


rule pirna_from_pirbase:
	input:
		bam="analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam",
		pirbase=config["pirbase"]
	output:
		"analyses/bowtie_mappings_genome_multi/txt_tables_per_file/" + version + "/{samplegroup}/" + x + "/{sample}." + x + ".txt" for x in ["pirbase"]
	shell:
		'''
            samplename=$(bin/sample_name.sh {input.bam})
			bin/rnaseq_tool.py {input.bam} -g {input.pirbase} --collapsed --filter -f piRNA -i name | awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' > {output}

		'''

rule mirna_and_precursor_from_mirBase:
	input:
		bam="analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam",
		mirbase=config["mirbase"]
	output:
		"analyses/bowtie_mappings_genome_multi/txt_tables_per_file/" + version + "/{samplegroup}/" + x + "/{sample}." + x + ".txt" for x in ["miRBasemiRNA","miRBaseprecursor"]
	shell:  
                '''
                samplename=$(bin/sample_name.sh {input.bam})

		bin/rnaseq_tool.py {input.bam} -g {input.mirbase} --collapsed --filter  -f miRNA -i Name |  awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' > {output[0]}
		bin/rnaseq_tool.py {input.bam} -g {input.mirbase} --collapsed --filter  -f miRNA_primary_transcript -i Name |  awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' > {output[1]}

                '''


rule mirna_and_precursor_from_mirgeneDB:
	input:
		bam="analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam",
		mirgenedb=config["mirgenedb"]
	
	output:
		"analyses/bowtie_mappings_genome_multi/txt_tables_per_file/" + version + "/{samplegroup}/" + x + "/{sample}." + x + ".txt" for x in ["mirgenedbprecursor","mirgenedbmirna"]

	shell:
		'''
        
		samplename=$(bin/sample_name.sh {input.bam})
		bin/rnaseq_tool.py {input.bam} -g {input.mirgenedb} --collapsed --filter -f miRNA_precursor -f miRNA_5p -f miRNA_3p -i ID |  awk -v s=$samplename -v precursor={output[0]} -v mirna={output[1]} 'BEGIN{{print"ID\t"s > precursor; print"ID\t"s > mirna}}{{if(/pre/) print >> precursor; else print >> mirna }}'

                '''


rule trna_from_gencode:
	input:
		bam="analyses/bowtie_mappings_genome_multi/{samplegroup}/{sample}.sorted.bam",
		trnagencode=config["trnagencode"]
	output:
		"analyses/bowtie_mappings_genome_multi/txt_tables_per_file/" + version + "/{samplegroup}/" + x + "/{sample}." + x + ".txt" for x in ["tRNAgencode"]
	
	shell:
		'''
        
		samplename=$(bin/sample_name.sh {input.bam})
		bin/rnaseq_tool.py {input.bam} -g {input.trnagencode} --collapsed --filter -f tRNA -i gene_id | awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' > {output}; #"gene_name to gene_id" 31.08.2020
		
        '''

rule create_final_count_tables:
	input:
		lambda wildcards : ["analyses/bowtie_mappings_genome_multi/txt_tables_per_file/" + wildcards.version + "/" + wildcards.samplegroup + "/" + wildcards.gene + "/" + x + "." + wildcards.gene + ".txt" for x in file_dict[wildcards.samplegroup]]

	output:
		"analyses/bowtie_mappings_genome_multi/count_tables/{version}/{samplegroup}/{gene}.count.table.csv"

	shell:
		'''
		bin/generic_table_creator.R {wildcards.gene} analyses/bowtie_mappings_genome_multi/txt_tables_per_file/{wildcards.version}/{wildcards.samplegroup}/{wildcards.gene}
		mv analyses/bowtie_mappings_genome_multi/txt_tables_per_file/{wildcards.version}/{wildcards.samplegroup}/{wildcards.gene}/{wildcards.gene}.count.table.csv analyses/bowtie_mappings_genome_multi/count_tables/{wildcards.version}/{wildcards.samplegroup}/{wildcards.gene}.count.table.csv
		'''



