import pandas as pd
configfile: "sunbeam_config.yml"
sample_file = config["all"]["samplelist_fp"]
samples = pd.read_csv(sample_file).iloc[:,1]
tags    = list(pd.read_csv(sample_file, header = None).iloc[:,0])

rule all:
    input:
        scaffolds=expand("amr/{tag}_metaspades/scaffolds.fasta", tag=tags),
        results=expand("amr/{tag}_staramr/results.xlsx", tag=tags),
        resistance=expand("amr/{tag}_collected_kraken_blast_results.tsv", tag=tags)

rule metaspades:
    input:
        r1="sunbeam_output/qc/cleaned/{tag}_1.fastq.gz",
        r2="sunbeam_output/qc/cleaned/{tag}_2.fastq.gz"
    output:
        metaspades=directory("amr/{tag}_metaspades"),
        scaffolds="amr/{tag}_metaspades/scaffolds.fasta"
    shell:
        """
        mkdir -p amr
        metaspades -t 10 -1 {input.r1} -2 {input.r2} -o {output.metaspades} -t 10
        """

rule staramr:
    input:
        metaspades = "amr/{tag}_metaspades/scaffolds.fasta"
    output:
        staramr = directory("amr/{tag}_staramr/"),
        results = "amr/{tag}_staramr/results.xlsx",
        resfinder = "amr/{tag}_staramr/resfinder.tsv"
    shell:
        """
        rm -rf {output.staramr}
        staramr search -o {output.staramr} {input.metaspades} -n 10
        """

rule select:
    input:
        metaspades = "amr/{tag}_metaspades/scaffolds.fasta",
        resistances = "amr/{tag}_staramr/resfinder.tsv"
    output:
        resistance_scaffolds = "amr/{tag}_resistance_scaffolds.fasta"
    shell:
        """
        python /media/deltagene/microbiome_2/AR_snakemake/scripts/select_scaffolds.fasta {input.metaspades} {input.resistances} {output.resistance_scaffolds}
        """

rule kraken:
    input:
        resistance_scaffolds = "amr/{tag}_resistance_scaffolds.fasta",
        kraken2_db = "/media/deltagene/microbiome_2/k2_plsupf_16gb_20240112"
    output:
        kraken2_out = "amr/{tag}_kraken2_results.txt"
    threads: 10
    shell:
        """
        kraken2 --output {output.kraken2_out} --db {input.kraken2_db} --threads {threads} {input.resistance_scaffolds}
        """
        
rule blast:
    input:
        resistance_scaffolds = "amr/{tag}_resistance_scaffolds.fasta",
        blast_db = "/media/deltagene/microbiome_2/k2_plsupf_16gb_20240112"
    output:
        blast_out = "amr/{tag}_blast_results.txt"
    threads: 10
    params:
        evalue = "1e-5"
    shell:
        """
        blastn -query {input.resistance_scaffolds} -db {input.blast_db} -out {output.blast_out} -evalue {params.evalue} -outfmt \"6 qseqid sseqid stitle pident length qstart qend sstart send evalue bitscore\" -num_threads {threads} -max_hsps 1 -max_target_seqs 5
        """
        
rule collect:
    input:
        resistance_scaffolds = "amr/{tag}_staramr/resfinder.tsv",
        kraken2_out = "amr/{tag}_kraken2_results.txt",
        blast_out = "amr/{tag}_blast_results.txt",
        taxonomy = "/media/data/Sunbeam/regi_kraken2/taxonomy"
    output:
        results = "amr/{tag}_collected_kraken_blast_results.tsv"
    shell:
        """
        python collect.py {input.resistance_scaffolds} {input.kraken2_out} {input.taxonomy} {input.blast_out} {output.results}
        """
"""
rule kraken_and_blast:
    input:
        scaffolds=expand("{tag}_metaspades/scaffolds.fasta",tag = tags),
        results=expand("{tag}_staramr/results.xlsx",tag = tags)
    output:
        "resistances_with_species_kraken_blast.tsv"
    shell:
        "python /media/deltagene/microbiome_2/AR_snakemake/scripts/collect_and_kraken.py ./"
"""
