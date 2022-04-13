'''
@author: Sinan U. Umu, PhD, sinanugur@gmail.com
'''

from collections import defaultdict
import os.path



rule all:
	input:
		expand("analyses/seqbuster/count_tables/{samplegroup}/{gene}.count.table.csv",samplegroup=set(directories),gene=gene),
		expand("analyses/seqbuster_new/{samplegroup}/result/{sample}.truncated_trimmed.mirna",zip,samplegroup=directories,sample=files),
		expand("analyses/mirtop/{samplegroup}/result/{sample}.truncated_trimmed.gff",zip,samplegroup=directories,sample=files)

rule seqcluster_collapse_or_unzip:
	input:
		"analyses/trimmed/{sample}.trimmed.fastq.gz"
	
	output:
		"analyses/seqbuster/collapse/{sample}.truncated_trimmed.fastq"

	params:
		out="analyses/seqbuster/collapse"
	
	run:
		shell("seqcluster collapse -f {input} -o {params.out}")


rule seqcluster_seqbuster:
	input:
		"analyses/seqbuster/collapse/{sample}.truncated_trimmed.fastq"
	
	output:
		"analyses/seqbuster/result/{sample}.truncated_trimmed.mirna"

	params:
		out="analyses/seqbuster/result",
		mirna="databases/miRNA.str",
		hairpin="databases/hsa-hairpin.fasta",
		sps="hsa"

	shell:
		"seqcluster seqbuster --out {params.out} --hairpin {params.hairpin} --mirna {params.mirna} --sps {params.sps} {input}"
rule seqcluster_seqbuster_new_output:
 	input:
 		"analyses/seqbuster/{samplegroup}/collapse/{sample}.truncated_trimmed.fastq"
 	output:
 		"analyses/seqbuster_new/{samplegroup}/result/{sample}.truncated_trimmed.mirna"

 	params:
 		prefix="analyses/seqbuster_new/{samplegroup}/result/{sample}.truncated_trimmed",
 		DB=mirbase_directory,
		parameters="-sub 1 -trim 3 -add 3 -s hsa"
 	shell:
 		"java -jar bin/miraligner.jar {params.parameters} -i {input} -db {params.DB} -o {params.prefix}"
		#"seqcluster seqbuster --out {params.out} --hairpin {params.hairpin} --mirna {params.mirna} --sps {params.sps} {input}"


rule mirtop:
	input:
		"analyses/seqbuster_new/{samplegroup}/result/{sample}.truncated_trimmed.mirna"
	output:
		"analyses/mirtop/{samplegroup}/result/{sample}.truncated_trimmed.gff"
	params:
		out="analyses/mirtop/{samplegroup}/result",
		prefix="analyses/seqbuster_new/{samplegroup}/result/{sample}.truncated_trimmed",
 		DB=mirbase_directory

	run:
		shell("""
		
		mirtop gff --format seqbuster --sps hsa --hairpin {params.DB}hairpin.fa  --gtf {params.DB}hsa.gff3 -o {params.out} {input}
		
		""")


rule create_seqbuster_count_tables:
	input:
		"analyses/seqbuster/{samplegroup}/result/{sample}.truncated_trimmed.mirna"
	
	output: 
		"analyses/seqbuster/{samplegroup}/result/{sample}.truncated_trimmed.mirna.hairpincounts.tsv",
		"analyses/seqbuster/{samplegroup}/result/{sample}.truncated_trimmed.mirna.mirnacounts.tsv",
		"analyses/seqbuster/{samplegroup}/result/{sample}.truncated_trimmed.mirna.isomircounts.tsv"
		

	shell:
		'''
		samplename=$(sample_name.sh {input});
		echo $samplename;
	        cat {input} | awk '!/name/{{print}}' | awk '{{print $1"\t"$3}}' | sort -k1,1 -u | awk -v s=$samplename 'BEGIN{{print"sequence\t"s}}{{print}}' > {output[2]};
        	cat {input} | awk '!/name/{{print}}' | awk -v  s=$samplename 'BEGIN{{print"miRNA\t"s}}{{a[$4]+=$3}}END{{for(i in a) print i"\t"a[i]}}' > {output[1]};
	        cat {input} | awk '!/name/{{print}}' | awk -v s=$samplename 'BEGIN{{print"hairpin\t"s}}{{a[$14]+=$3}}END{{for(i in a) print i"\t"a[i]}}' > {output[0]};

		'''
rule create_final_seqbuster_tables:
	input:
                lambda wildcards : ["analyses/seqbuster/" + wildcards.samplegroup + "/result/" + x + ".truncated_trimmed.mirna." + wildcards.gene + "counts" + ".tsv" for x in file_dict[wildcards.samplegroup]]

	output:
		"analyses/seqbuster/count_tables/{samplegroup}/{gene}.count.table.csv"

	shell:
		'''
		generic_table_creator.R {wildcards.gene}counts analyses/seqbuster/{wildcards.samplegroup}/result/
                
		mv analyses/seqbuster/{wildcards.samplegroup}/result/{wildcards.gene}counts.count.table.csv analyses/seqbuster/count_tables/{wildcards.samplegroup}/{wildcards.gene}.count.table.csv


		'''






		

 
