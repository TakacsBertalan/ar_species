import pandas as pd
import os

#Default path to configfile. Needs to be modified before running the pipeline or a new configfile needs to be set with the -c flag
configfile: config.get("configfile", "/media/deltagene/microbiome_2/AR_snakemake/config.yaml")
print(config)

sample_file = config["all"]["samplelist_fp"]
output_folder = config["all"]["output_folder"]
print(sample_file)

tags = list(pd.read_csv(sample_file, header = None).iloc[:,0])

try:
	os.mkdir({output_folder})
except:
	print(f"Output folder {output_folder} already exists")

rule all:
    input:
        scaffolds=expand(output_folder + "/{tag}_metaspades/scaffolds.fasta", tag=tags),
        results=expand(output_folder + "/{tag}_staramr/results.xlsx", tag=tags),
        resistance=expand(output_folder + "/{tag}_collected_kraken_blast_results.tsv", tag=tags)

rule metaspades:
    input:
        r1="{config[all][input_fp]}/{tag}_1.fastq.gz",
        r2="{config[all][input_fp]}/{tag}_2.fastq.gz"
    output:
        metaspades=directory("{output_folder}/{tag}_metaspades"),
        scaffolds="{output_folder}/{tag}_metaspades/scaffolds.fasta"
    threads:
        config['metaspades']['threads']
    shell:
        """
        mkdir -p amr
        metaspades -t 10 -1 {input.r1} -2 {input.r2} -o {output.metaspades} -t {threads}
        """

rule staramr:
    input:
        metaspades = "amr/{tag}_metaspades/scaffolds.fasta"
    output:
        staramr = directory("{output_folder}/{tag}_staramr/"),
        results = "{output_folder}/{tag}_staramr/results.xlsx",
        resfinder = "{output_folder}/{tag}_staramr/resfinder.tsv"
    threads:
        config['staramr']['threads']
    shell:
        """
        rm -rf {output.staramr}
        staramr search -o {output.staramr} {input.metaspades} -n {threads}
        """

rule select:
    input:
        metaspades = "{output_folder}/{tag}_metaspades/scaffolds.fasta",
        resistances = "{output_folder}/{tag}_staramr/resfinder.tsv"
    output:
        resistance_scaffolds = "{output_folder}/{tag}_resistance_scaffolds.fasta"
    shell:
        """
        python {config[all][scripts]}/select_scaffolds.py {input.metaspades} {input.resistances} {output.resistance_scaffolds}
        """

rule kraken:
    input:
        resistance_scaffolds = "{output_folder}/{tag}_resistance_scaffolds.fasta",
        kraken2_db = config["kraken2"]["db"]
    output:
        kraken2_out = "{output_folder}/{tag}_kraken2_results.txt"
    threads:
        config['kraken2']['threads']
    shell:
        """
        kraken2 --output {output.kraken2_out} --db {input.kraken2_db} --threads {threads} {input.resistance_scaffolds}
        """
        
rule blast:
    input:
        resistance_scaffolds = "{output_folder}/{tag}_resistance_scaffolds.fasta",
        
    output:
        blast_out = "{output_folder}/{tag}_blast_results.txt"
    threads:
        config['blast']['threads']
    params:
        evalue = config["blast"]["e-value"],
        blast_db = config["blast"]["db"]
    shell:
        """
        blastn -query {input.resistance_scaffolds} -db {params.blast_db} -out {output.blast_out} -evalue {params.evalue} -outfmt \"6 qseqid sseqid stitle pident length qstart qend sstart send evalue bitscore\" -num_threads {threads} -max_hsps 1 -max_target_seqs 5
        """
        
rule collect:
    input:
        resistance_scaffolds = "{output_folder}/{tag}_staramr/resfinder.tsv",
        kraken2_out = "{output_folder}/{tag}_kraken2_results.txt",
        blast_out = "{output_folder}/{tag}_blast_results.txt",
        taxonomy = "/media/data/Sunbeam/regi_kraken2/taxonomy"
    output:
        results = "{output_folder}/{tag}_collected_kraken_blast_results.tsv"
    shell:
        """
        python {config[all][scripts]}/improved_collect.py {input.resistance_scaffolds} {input.kraken2_out} {input.taxonomy} {input.blast_out} {output.results}
        """
