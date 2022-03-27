from collections import defaultdict
#directories, files, = glob_wildcards("data/raw/{samplegroup}/{sample}.fastq.gz")

#files, = glob_wildcards("data/{sample}.fastq.gz")

include: "bowtie2.smk"
#gene=['tRF']


rule trna_derived_by_MINTmap:
        input: "analyses/trimmed/{sample}.trimmed.fastq.gz"

        output: "analyses/MINTmap/{sample}-MINTmap_v1-exclusive-tRFs.expression.txt",
                "analyses/MINTmap/{sample}-MINTmap_v1-ambiguous-tRFs.expression.txt",
		"analyses/MINTmap/txt_tables_per_file/{sample}.MINTmap.txt"
        params:
                p="analyses/MINTmap/{sample}",
                l="databases/LookupTable.tRFs.MINTmap_v1.txt",
                s="databases/tRNAspace.Spliced.Sequences.MINTmap_v1.fa",
                o="databases/OtherAnnotations.MINTmap_v1.txt",
                j="databases/MINTplates/"


        shell:
                '''
		input_name=$(basename {input})
		samplename=${{input_name/.trimmed.fastq.gz/}}

                MINTmap.pl -p {params.p} -f {input} -l {params.l} -s {params.s} -o {params.o} -j {params.j};

		cat {output[0]} {output[1]} | gawk -v s=$samplename 'BEGIN{{print "ID\t"s}}!/MINTbase/{{print $1"\t"$4}}' > {output[2]};

                '''


rule create_tRF_count_tables:
        input:
                lambda wildcards : ["analyses/MINTmap/txt_tables_per_file/" + s + ".MINTmap.txt" for s in files]

        output:
                "results/count_tables/tRF.tsv"

        shell:
                '''
                #Rscript --no-environ workflow/scripts/generic_table_creator.R MINTmap analyses/MINTmap/txt_tables_per_file/
                #mv analyses/MINTmap/txt_tables_per_file/MINTmap.tsv results/count_tables/tRF.tsv

                echo "{output} {input}" | xargs Rscript --no-environ ./workflow/scripts/sncrna_table_creator.R

                '''


