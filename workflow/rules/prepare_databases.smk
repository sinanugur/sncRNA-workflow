from collections import defaultdict





rule prepare_gencode:
    input:
        "databases/gencode.gff3.gz"
    output:
        temp("databases/gencode.gff3")
    shell:
        "gunzip -c {input} > gencode.gff3"

rule prepare_pirna:
    input:
        "databases/piRNA.gtf.gz"
    output:
        temp("databases/piRNA.gtf")
    shell:
        "gunzip -c {input} > piRNA.gtf"