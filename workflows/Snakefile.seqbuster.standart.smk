from collections import defaultdict
import os.path


if config.get('directory') is None:
	directories, files, = glob_wildcards("data/raw/{samplegroup}/{sample}.fastq.gz")
	#glob_wildcards("WGS/tests/fastq/{group}/{sample}.fastq.gzadr_R1.fastq")
else: #a directory name can be provided for a single group run
	files, = glob_wildcards("data/raw/" + config.get('directory') + "/{sample}.fastq.gz")
	directories=[config.get('directory')]*len(files)



mirbase_directory="/WORKING/databases/miRBase/"




gene=["hairpin","isomir","mirna"]


file_dict=defaultdict(list)

#assign files to correct directories
for i in zip(directories,files):
        file_dict[i[0]].append(i[1])


print(set(directories))


rule all:
	input:
		expand("analyses/seqbuster/count_tables/{samplegroup}/{gene}.count.table.csv",samplegroup=set(directories),gene=gene),
		expand("analyses/seqbuster_new/{samplegroup}/result/{sample}.truncated_trimmed.mirna",zip,samplegroup=directories,sample=files),
		expand("analyses/mirtop/{samplegroup}/result/{sample}.truncated_trimmed.gff",zip,samplegroup=directories,sample=files)

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




rule seqcluster_collapse_or_unzip:
	input:
		"data/trimmed/{samplegroup}/{sample}.truncated.fastq.gz"
	
	output:
		"analyses/seqbuster/{samplegroup}/collapse/{sample}.truncated_trimmed.fastq"

	params:
		o="analyses/seqbuster/{samplegroup}/collapse"
	
	run:
		#shell("/WORKING/apps/miniconda2/envs/python2env/bin/seqcluster collapse -f {input} -o {params.o}")
		shell("seqcluster collapse -f {input} -o {params.o}")


rule seqcluster_seqbuster:
	input:
		"analyses/seqbuster/{samplegroup}/collapse/{sample}.truncated_trimmed.fastq"
	
	output:
		"analyses/seqbuster/{samplegroup}/result/{sample}.truncated_trimmed.mirna"

	params:
		out="analyses/seqbuster/{samplegroup}/result",
		mirna=mirbase_directory + "miRNA.str",
		hairpin=mirbase_directory + "hsa-hairpin.fasta",
		sps="hsa"

	shell:
		#"/WORKING/apps/miniconda2/envs/python2env/bin/seqcluster seqbuster --out {params.out} --hairpin {params.hairpin} --mirna {params.mirna} --sps {params.sps} {input}"
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






		

 
