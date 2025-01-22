import os
import subprocess
import csv
from sys import argv

def collect_scaffolds(resfinder, metaspades, res_scaffolds)
    entries = {}
	scaffolds = {}
	i = 0
    with open(resfinder, "r") as res:
		for line in res:
			if first_line:
				first_line = False
			else:
				comp = line.split("\t")
				if "unknown" not in comp[2]:
				entries[str(i)] = [name, comp[1], comp[2], comp[3], comp[4], comp[6], comp[10].rstrip()]
				scaffolds[comp[6]] = ""
				i += 1
	with open(metaspades, "r") as scaf:
	    read_in = False
		read_string = ""
		found = ""
		for line in scaf:
	    	if line[0] == ">":
				read_in = False
				if found != "":
	    			scaffolds[found] = read_string
		    		found = ""
					read_string = ""
			    	if line.rstrip()[1:] in list(scaffolds.keys()):
						found = line.rstrip()[1:]
						read_in = True
					elif read_in:
						read_string += line.rstrip()
				
				if found != "":
					scaffolds[found] = read_string
					
	with open(res_scaffolds, "w") as fasta:
		for key in entries:
			if "unknown" not in entries[key][2]:
				fasta_string = ">" + entries[key][0] + "_" + entries[key][1] + "_" + entries[key][5]
				fasta.write(fasta_string + "\n")
				fasta.write(scaffolds[entries[key][-2]] + "\n")

collect_scaffolds(argv[1],argv[2],argv[3])
