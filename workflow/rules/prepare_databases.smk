from collections import defaultdict

import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()


rule a:
	input:
		multiext(
            "databases/hg38",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        )

rule prepare_gencode:
    input:
        "databases/gencode.gff3.gz"
    output:
        temp("databases/gencode.gff3")
    shell:
        "gunzip -c {input} > databases/gencode.gff3"

rule prepare_pirna:
    input:
        "databases/piRNA.gtf.gz"
    output:
        temp("databases/piRNA.gtf")
    shell:
        "gunzip -c {input} > databases/piRNA.gtf"


rule prepare_bowtie2_index:
    input:
        "databases/hg38.fa.gz"	
    output:
        multiext(
            "databases/hg38",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
        
    conda:
        "../envs/main.yaml"
    threads: 15

    shell:
        "bowtie2-build {input} databases/hg38"
