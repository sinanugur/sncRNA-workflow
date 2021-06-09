import os

directories, files, = glob_wildcards("data/raw/{samplegroup}/{sample}.fastq.gz")

configfile: "workflows/config.yaml"
version=str(config["version"])

directory_bowtie2="analyses/bowtie_mappings_genome_multi/count_tables/" + version + "/"
directory_seqbuster="analyses/seqbuster/count_tables/"
directory_results="results/" + version + "/"
directory_trf="analyses/MINTmap/count_tables/"



print("Current working directory: " + os.getcwd())
directory_workdir=os.getcwd()

database_vector=["gene.gencode","miRBasemiRNA","miRBaseprecursor","mirgenedbmirna","mirgenedbprecursor","pirbase","tRNAgencode","exons.gencode","gencode.uniq.exon","pirbase.uniq","tRNAgencode.uniq"]

rule all:
	input:
		#expand([directory_results + "{samplegroup}/" + x + ".count.table.csv" for x in database_vector],samplegroup=set(directories)),
		#expand([directory_results + "{samplegroup}/" + x + ".count.table.csv" for x in ["mirna","hairpin","isomir","MINTmap"]],samplegroup=set(directories)),
		[directory_results  + x + ".combined.table.csv" for x in database_vector],
		[directory_results  + x + ".combined.table.csv" for x in ["mirna","isomir","hairpin","MINTmap"]],
		[directory_results + x + ".combined.table.csv" for x in ["lincRNA","misc_RNA","protein_coding","snoRNA","snRNA","scaRNA","antisense"]],
		directory_results  + "mirjanus_csv_v" + version + ".tar.gz"
		



rule copy_tables_to_results:
	input:
		[directory_bowtie2 + "{samplegroup}/" + x + ".count.table.csv" for x in database_vector],
		directory_trf + "{samplegroup}/MINTmap.count.table.csv",
		[directory_seqbuster + "{samplegroup}/" + x for x in ["mirna.count.table.csv","isomir.count.table.csv","hairpin.count.table.csv"]]

	output:
		[directory_results + "{samplegroup}/" + x + ".count.table.csv" for x in database_vector],
		directory_results + "{samplegroup}/MINTmap.count.table.csv",
		[directory_results + "{samplegroup}/" + x for x in ["mirna.count.table.csv","isomir.count.table.csv","hairpin.count.table.csv"]]

	run:
		for i in input:
			shell("cp {i} {directory_results}/{wildcards.samplegroup}")

rule create_tables:
	input:
		expand([directory_results + "{samplegroup}/" + x + ".count.table.csv" for x in database_vector],samplegroup=set(directories)),
		expand(directory_results + "{samplegroup}/MINTmap.count.table.csv",samplegroup=set(directories)),
		expand([directory_results + "{samplegroup}/" + x for x in ["mirna.count.table.csv","isomir.count.table.csv","hairpin.count.table.csv"]],samplegroup=set(directories))
	output:
		[directory_results  + x + ".combined.table.csv" for x in database_vector],
		[directory_results  + x + ".combined.table.csv" for x in ["mirna","isomir","hairpin","MINTmap"]]
	run:
		shell("""
		cd {directory_results}
		{directory_workdir}/bin/build_final_table_databases.R

		""")


rule separate_gene_types:
	input:
		directory_results + "exons.gencode.combined.table.csv"
	output:
		[directory_results + x + ".combined.table.csv" for x in ["lincRNA","misc_RNA","protein_coding","snoRNA","snRNA","scaRNA","antisense"]]
	run:
		shell("""
		cd {directory_results}
		{directory_workdir}/bin/separate_gene_types.R

		""")


rule create_csv_zip:
	input:
		expand([directory_results + "{samplegroup}/" + x + ".count.table.csv" for x in database_vector],samplegroup=set(directories)),
		expand(directory_results + "{samplegroup}/MINTmap.count.table.csv",samplegroup=set(directories)),
		expand([directory_results + "{samplegroup}/" + x for x in ["mirna.count.table.csv","isomir.count.table.csv","hairpin.count.table.csv"]],samplegroup=set(directories)),
		[directory_results + x + ".combined.table.csv" for x in ["lincRNA","misc_RNA","protein_coding","snoRNA","snRNA","scaRNA","antisense"]]
		
	output:
		directory_results  + "mirjanus_csv_v" + version + ".tar.gz"
	run:
		shell("""
		cd {directory_results}
		find -name "*csv" | tar -czvf `basename {output}` -T -
		""")