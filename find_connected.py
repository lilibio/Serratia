import sys

def extract_genome_gene_id(gene):
    segments = gene.split('_')
    genome_id = '_'.join(segments[:-1])
    gene_number = segments[-1]
    return genome_id, gene_number

def find_connected_genes(gene_list):
    connected_groups = []
    current_group = []

    for gene in sorted(gene_list):
        if not current_group or int(extract_genome_gene_id(gene)[1]) == int(extract_genome_gene_id(current_group[-1])[1]) + 1:
            current_group.append(gene)
        else:
            if len(current_group) > 1:  # Only include groups with more than one gene
                connected_groups.append(current_group)
            current_group = [gene]

    if current_group and len(current_group) > 1:
        connected_groups.append(current_group)

    return connected_groups

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py gene_list_file.txt")
        sys.exit(1)

    gene_list_file = sys.argv[1]

    with open(gene_list_file, 'r') as file:
        genes = [line.strip() for line in file]

    connected_groups = find_connected_genes(genes)

    for group in connected_groups:
        print("Connected Genes:", ', '.join(group))
