from collections import defaultdict
directories, files, = glob_wildcards("data/raw/{samplegroup}/{sample}.fastq.gz")

annotation=["gencode"]


gene=['MINTmap']


file_dict=defaultdict(list)

#assign files to correct directories
for i in zip(directories,files):
        file_dict[i[0]].append(i[1])



rule all:
	input:
		expand("analyses/MINTmap/count_tables/{samplegroup}/{gene}.count.table.csv",samplegroup=set(directories),gene=gene)


rule adapter_trimming:
	input:
		"data/raw/{samplegroup}/{sample}.fastq.gz"

	output:
		"data/trimmed/{samplegroup}/{sample}.truncated.fastq.gz"
	
	threads: 15

	shell:
		# "Adapter_removal.sh {input} {output}"
		"""
		input={input}
		basename=${input/.fastq.gz/}
		AdapterRemoval --file1 {input} --basename $basename --trimns --trimqualities --gzip --threads {threads} --minlength 17 --minquality 30;
		mv "$basename".truncated.gz "$basename".truncated.fastq.gz;
		mv "$basename".discarded.gz "$basename".discarded.fastq.gz;
		"""


rule trna_derived_by_MINTmap:
        input: "data/trimmed/{samplegroup}/{sample}.truncated.fastq.gz"

        output: "analyses/MINTmap/{samplegroup}/{sample}-MINTmap_v1-exclusive-tRFs.expression.txt",
                "analyses/MINTmap/{samplegroup}/{sample}-MINTmap_v1-ambiguous-tRFs.expression.txt",
		"analyses/MINTmap/txt_tables_per_file/{samplegroup}/{sample}.MINTmap.txt"
        params:
                p="analyses/MINTmap/{samplegroup}/{sample}",
                l="/WORKING/databases/mintmap_tRF/MINTmap-master/LookupTable.tRFs.MINTmap_v1.txt",
                s="/WORKING/databases/mintmap_tRF/MINTmap-master/tRNAspace.Spliced.Sequences.MINTmap_v1.fa",
                o="/WORKING/databases/mintmap_tRF/MINTmap-master/OtherAnnotations.MINTmap_v1.txt",
                j="/WORKING/databases/mintmap_tRF/MINTmap-master/MINTplates/"


        shell:
                '''
		samplename=$(sample_name.sh {input});

                MINTmap.pl -p {params.p} -f {input} -l {params.l} -s {params.s} -o {params.o} -j {params.j};

		cat {output[0]} {output[1]} | gawk -v s=$samplename 'BEGIN{{print "ID\t"s}}!/MINTbase/{{print $1"\t"$4}}' > {output[2]};

                '''


rule create_final_count_tables:
        input:
                lambda wildcards : ["analyses/MINTmap/txt_tables_per_file/" + wildcards.samplegroup + "/" + x + "." + wildcards.gene + ".txt" for x in file_dict[wildcards.samplegroup]]

        output:
                "analyses/MINTmap/count_tables/{samplegroup}/{gene}.count.table.csv"

        shell:
                '''
                generic_table_creator.R {wildcards.gene} analyses/MINTmap/txt_tables_per_file/{wildcards.samplegroup}/
                mv analyses/MINTmap/txt_tables_per_file/{wildcards.samplegroup}/{wildcards.gene}.count.table.csv analyses/MINTmap/count_tables/{wildcards.samplegroup}/{wildcards.gene}.count.table.csv
                '''


