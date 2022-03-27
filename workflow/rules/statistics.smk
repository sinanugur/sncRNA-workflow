from collections import defaultdict

'''
@author: Sinan U. Umu, PhD, sinanugur@gmail.com
'''
include: "bowtie2.smk"
include: "prepare_databases.smk"

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4


files, = glob_wildcards("data/{sample}.fastq.gz")

configfile: "config/config.yaml"

rule uniq_mapping_gencode:
        input:
                bam="analyses/bowtie_mappings_genome_multi/{sample}.sorted.bam",
                gencode="databases/gencode.gff3"
        output:
                "analyses/statistics/txt_tables_per_file/" + x + "/{sample}." + x + ".txt" for x in ["gencode.uniq.exon"]
        params:
                par="--collapsed --unique -f exon -i gene_type"
        shell:
                """

                                input_name=$(basename {input.bam})
                samplename=${{input_name/.sorted.bam/}}
                workflow/scripts/rnaseq_tool.py  {input.bam} -g {input.gencode} {params.par} | awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' >  {output}
                """

rule uniq_mapping_trna:
        input:
                bam="analyses/bowtie_mappings_genome_multi/{sample}.sorted.bam",
                trnagencode="databases/tRNA.forstatistics.gff3"

        output:
                "analyses/statistics/txt_tables_per_file/" + x + "/{sample}." + x + ".txt" for x in ["tRNA.uniq"]
        params:
                par = "--collapsed --unique -f tRNA -i stat_type"

        shell:
                """
                                input_name=$(basename {input.bam})
                samplename=${{input_name/.sorted.bam/}}
                workflow/scripts/rnaseq_tool.py {input.bam} -g {input.trnagencode} {params.par} | awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' >  {output}
                """




rule uniq_mapping_pirna:
        input:
                bam="analyses/bowtie_mappings_genome_multi/{sample}.sorted.bam",
                pirbase="databases/piRNA.gtf"
        output:
                "analyses/statistics/txt_tables_per_file/" + x + "/{sample}." + x + ".txt" for x in ["piRNA.uniq"]
        params:
                par = "--collapsed --unique -f piRNA -i gene_type"

        shell:
                """
                                input_name=$(basename {input.bam})
                samplename=${{input_name/.sorted.bam/}}
                workflow/scripts/rnaseq_tool.py  {input.bam} -g {input.pirbase} {params.par} | awk -v s=$samplename 'BEGIN{{print"ID\t"s}}{{print}}' >  {output}
        """


rule create_stats_count_tables:
        input:
                lambda wildcards : ["analyses/statistics/txt_tables_per_file/" + wildcards.gene + "/" + s + "." + wildcards.gene + ".txt" for s in files]

        output:
                "results/statistics/{gene}.tsv"

        conda:
                "../envs/main.yaml"

        shell:
                '''
                echo "{output} {input}" | xargs Rscript --no-environ ./workflow/scripts/sncrna_table_creator.R

                
                '''
