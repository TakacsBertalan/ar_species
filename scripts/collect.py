from sys import argv

def translate_kraken(kraken2_output, kraken2_taxonomy):
	
	print("Collecting Kraken2 results")

	# Paths to files
	taxonomy_nodes = kraken2_taxonomy + "/nodes.dmp"
	taxonomy_names = kraken2_taxonomy + "/names.dmp"

	# Step 1: Load taxonomy names and ranks into dictionaries
	tax_id_to_name = {}
	tax_id_to_parent = {}
	tax_id_to_rank = {}

	tax_dict = {}

	# Load `names.dmp`
	with open(taxonomy_names, "r") as names_file:
		for line in names_file:
			parts = line.strip().split("\t|\t")
			tax_id = parts[0].strip()
			name = parts[1].strip()
			tax_id_to_name[tax_id] = name
			
	# Load `nodes.dmp`
	with open(taxonomy_nodes, "r") as nodes_file:
		for line in nodes_file:
			parts = line.strip().split("\t|\t")
			tax_id = parts[0].strip()
			parent_id = parts[1].strip()
			rank = parts[2].strip()
			tax_id_to_parent[tax_id] = parent_id
			tax_id_to_rank[tax_id] = rank
	# Step 2: Helper function to construct lineage up to species or higher

	def get_full_lineage(tax_id):
		lineage = []
		current_id = tax_id
		seen_ids = set()  # Prevent infinite loop

		while current_id in tax_id_to_parent:
			# Prevent infinite loop
			if current_id in seen_ids:
				break
			seen_ids.add(current_id)

			rank = tax_id_to_rank.get(current_id, "Unknown")
			name = tax_id_to_name.get(current_id, "Unknown Taxon")

			# Stop at root or if no parent is found
			if current_id == "1" or rank == "no rank":
				lineage .append(name)
				break

			# Only add ranks we're interested in
			if rank in ["species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"]:
				lineage.append(name)

			# Move to parent
			current_id = tax_id_to_parent[current_id]

		return " > ".join(reversed(lineage)) if lineage else "Unknown Taxon"

	# Step 3: Process Kraken2 output and filter by species level or higher
	with open(kraken2_output, "r") as infile:
		for line in infile:
			parts = line.strip().split("\t", maxsplit=4)
			classification = "Classified" if parts[0] == "C" else "Unclassified"
			seq_id = parts[1]
			tax_id = parts[2]
			tax_dict[parts[1]] = [parts[2:]]

			if classification == "Classified":
				tax_dict[parts[1]] = ["Classified"]

				rank = tax_id_to_rank.get(tax_id, "Unknown")
				full_taxonomic_name = get_full_lineage(tax_id)

				# Include only species level or higher
				if rank in ["species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"]:
					tax_dict[seq_id] = [seq_id, classification, tax_id, full_taxonomic_name, rank]
			else:
				tax_dict[seq_id] = [seq_id, classification, "-", "Unclassified", "-"]

	return tax_dict

def read_blast_result(blast_outp):
	results = {}
	with open(blast_outp, "r") as blast:
		for line in blast:

			if line.split("\t")[0] not in results:
				results[line.split("\t")[0]] = [line.split("\t")[2]]
			else:
				results[line.split("\t")[0]].append(line.split("\t")[2])
	return results

def collect_results(resistances, kraken2_output, kraken2_taxonomy, blast_output, collected_results):
	entries = {}
	res_dict = {}
	scaffolds = {}
	i = 0
	first_line = True
	name = resistances.split("/")[-2].split("_")[0]
	with open(resistances, "r") as res:
		for line in res:
			if first_line:
				first_line = False
			else:
				comp = line.split("\t")
				if "unknown" not in comp[2]:
					entries[str(i)] = [name, comp[1], comp[2], comp[3], comp[4], comp[6], comp[10].rstrip()]
					i += 1
	
	with open(collected_results, "w") as collection:
		collection.write("Sample\tGene\tPredicted phenotype\t%Identity\t%Overlap\tScaffold name\tKraken2 Classified?\tTax ID\tLineage\tLowest identified species\tBest 5 BLAST hits\n")
		kraken_results = translate_kraken(kraken2_output, kraken2_taxonomy)
		blast_results = read_blast_result(blast_output)
		for key in entries:
			collection.write("\t".join(entries[key][0:5]) + "\t")

			if len(kraken_results[entries[key][0] + "_" + entries[key][1] + "_" + entries[key][5]]) == 1:
				collection.write(entries[key][0] + "_" + entries[key][1] + "_" + entries[key][5] + "\tna\tna\tna\tna\t")
			else:
				collection.write("\t".join(kraken_results[entries[key][0] + "_" + entries[key][1] + "_" + entries[key][5]]) + "\t")

			collection.write(", ".join(blast_results[entries[key][0] + "_" + entries[key][1] + "_" + entries[key][5]]))
			collection.write("\n")
			
collect_results(argv[1], argv[2], argv[3], argv[4], argv[5])
