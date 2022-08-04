from collections import defaultdict
#directories, files, = glob_wildcards("data/raw/{samplegroup}/{sample}.fastq.gz")

#files, = glob_wildcards("data/{sample}.fastq.gz")

include: "bowtie2.smk"
#gene=['tRF']

'''
@author: Sinan U. Umu, PhD, sinanugur@gmail.com
'''


rule trna_derived_by_MINTmap:
        input: "analyses/trimmed/{sample}.trimmed.fastq.gz"

        output: "analyses/MINTmap/{sample}-MINTmap_v1-exclusive-tRFs.expression.txt",
                "analyses/MINTmap/{sample}-MINTmap_v1-ambiguous-tRFs.expression.txt"
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
                '''


rule trf_txt_tables:
    input:
        "analyses/MINTmap/{sample}-MINTmap_v1-exclusive-tRFs.expression.txt",
        "analyses/MINTmap/{sample}-MINTmap_v1-ambiguous-tRFs.expression.txt"
    output:
        "analyses/MINTmap/txt_tables_per_file/{sample}.MINTmap.txt"
    shell:
        """
        input_name=$(basename {input[0]})
	samplename=${{input_name/-MINTmap_v1-exclusive-tRFs.expression.txt/}}

        cat {input[0]} {input[1]} | gawk -v s=$samplename 'BEGIN{{print "ID\t"s}}!/MINTbase/{{print $1"\t"$4}}' > {output};



        """


rule create_tRF_count_tables:
        input:
                lambda wildcards : ["analyses/MINTmap/txt_tables_per_file/" + s + ".MINTmap.txt" for s in files]

        output:
                "results/count_tables/tRF.tsv"

        shell:
                '''

		for i in {input}; do echo $i; done > analyses/MINTmap/tmp.csv
		Rscript --no-environ ./workflow/scripts/sncrna_table_creator.R {output} analyses/MINTmap/tmp.csv

                rm analyses/MINTmap/tmp.csv

                '''


