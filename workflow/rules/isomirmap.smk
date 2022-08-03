from collections import defaultdict
#directories, files, = glob_wildcards("data/raw/{samplegroup}/{sample}.fastq.gz")

#files, = glob_wildcards("data/{sample}.fastq.gz")

include: "bowtie2.smk"
#gene=['tRF']

'''
@author: Sinan U. Umu, PhD, sinanugur@gmail.com
'''


rule isomir_derived_by_isomirmap:
        input:
            "analyses/trimmed/{sample}.trimmed.fastq.gz"

        output: 
            "analyses/isoMiRmap/{sample}-IsoMiRmap_v5-exclusive-isomiRs.expression.txt",
            "analyses/isoMiRmap/{sample}-IsoMiRmap_v5-ambiguous-isomiRs.expression.txt",
            "analyses/isoMiRmap/{sample}-IsoMiRmap_v5-snps-isomiRs.expression.txt",
            "analyses/isoMiRmap/txt_tables_per_file/{sample}.isomirmap.txt"

        shell:
            """
		input_name=$(basename {input})
		samplename=${{input_name/.trimmed.fastq.gz/}}

        isoMiRmap/IsoMiRmap.py {input} --m ./isoMiRmap/MappingBundles/miRBase/ --p analyses/isoMiRmap/"$samplename"


		cat {output[0]} {output[1]} | grep -v "^#" | gawk -v s=$samplename 'BEGIN{{print "ID\t"s}}!/License/{{print $1"\t"$4}}' > {output[3]};

        cat {output[2]} | grep -v "^#" | gawk '!/License/{{print $1"\t"$3}}' >> {output[3]};

            """


rule create_isomirmap_count_tables:
        input:
            lambda wildcards : ["analyses/isoMiRmap/txt_tables_per_file/" + s + ".isomirmap.txt" for s in files]

        output:
            "results/count_tables/isomiR.tsv"

        shell:
            """
		files=$(for i in {input}; do printf $i","; done)
                
		Rscript --no-environ ./workflow/scripts/sncrna_table_creator.R {output} $files
            """


