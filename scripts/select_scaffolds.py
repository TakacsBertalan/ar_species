import os
import subprocess
import csv
from sys import argv
from Bio import SeqIO

def collect_scaffolds(metaspades, resfinder, res_scaffolds):
	entries = {}
	scaffolds = {}
	i = 0
	first_line = True
	name = metaspades.split("/")[-2].split("_")[0]
	with open(resfinder, "r") as res:
		for line in res:
			if first_line:
				first_line = False
			else:

				comp = line.split("\t")
				print(comp)
				if "unknown" not in comp[2]:
					entries[str(i)] = [name, comp[1], comp[2], comp[3], comp[4], comp[6], comp[10].rstrip()]
					scaffolds[comp[6]] = ""
				i += 1

	for record in SeqIO.parse(metaspades, "fasta"):
		if record.id in scaffolds:
			scaffolds[record.id] = str(record.seq)
					
	with open(res_scaffolds, "w") as fasta:
		for key in entries:
			if "unknown" not in entries[key][2]:
				fasta_string = ">" + entries[key][0] + "_" + entries[key][1] + "_" + entries[key][5]
				fasta.write(fasta_string + "\n")
				fasta.write(scaffolds[entries[key][-2]] + "\n")

print("collect_scaffolds " + "\t".join([argv[1],argv[2],argv[3]]))
collect_scaffolds(argv[1],argv[2],argv[3])
