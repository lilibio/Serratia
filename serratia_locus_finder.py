import sys

def read_input_file(file_path):
    data = []
    with open(file_path, 'r') as file:
        lines = file.readlines()[1:]
        for line in lines:
            columns = line.strip().split('\t')
            
            # Extract genes from the fifth column separated by commas
            genes_fifth_column = columns[3].split(',')
            
            # Iterate over each gene in the fifth column
            for gene_fifth_column in genes_fifth_column:
                # Extract the second part after the underscore of the first gene in the fourth column
                first_gene_parts = columns[4].split(',')
                if first_gene_parts:
                    first_gene_part = first_gene_parts[0].split('_')[0]
                    if len(first_gene_part) > 1:
                      data.append((first_gene_part, gene_fifth_column))
    
    return data


def is_connected(gene1, gene2):
    gene_number1 = int(gene1.split('_')[-1][1:])
    gene_number2 = int(gene2.split('_')[-1][1:])
    distance = abs(gene_number1 - gene_number2)
    return distance <= 10

def get_subject_gene(query_gene, data):
    for subject_gene, q_gene in data:
        if q_gene == query_gene:
            return subject_gene
    return None

def count_core_components_in_column(column_data):
    core_components = [
        'TssA', 'TssB', 'TssC', 'TssD', 'TssE', 'TssF', 'TssG', 'TssH', 'TssI','TssK', 'TssL', 'TssM','TssJ'
    ]
    present_core_components = [component for gene in column_data for component in core_components if component in gene]
    core_count = len(set(present_core_components))
    return core_count, present_core_components

def sort_genes_by_number(fifth_column_data):
    sorted_genes = sorted(fifth_column_data, key=lambda gene: int(gene.split('_')[-1][1:]))
    return sorted_genes

def extract_first_gene_part(column_data):
    # Split the first gene in the fourth column by commas
    genes = column_data[0].split(',')
    # Extract the first part of the first gene before the underscore
    if genes:
        first_gene = genes[0].split('_')[0]
        return first_gene
    else:
        return None


def main():
    if len(sys.argv) != 2:
        return
    input_file_path = sys.argv[1]
    data = read_input_file(input_file_path)
    fourth_column_data = [item[0] for item in data]
    fifth_column_data = [item[1] for item in data]   
    core_count, present_core_components = count_core_components_in_column(fourth_column_data)
    sorted_fifth_column_data = sort_genes_by_number(fifth_column_data)
    print(f"Number of core components: {core_count}")
    print(f"{present_core_components}")

    loci = []
    current_locus = []   
    for gene in sorted_fifth_column_data:
        if not current_locus:
            current_locus.append(gene)
        else:
            last_gene = current_locus[-1]
            if is_connected(last_gene, gene):
                current_locus.append(gene)
            else:
                if len(current_locus) >= 3:
                    loci.append(current_locus)
                current_locus = [gene]    
    if len(current_locus) >= 3:
        loci.append(current_locus)
    print("Number of T6SS loci:", len(loci))

    for i, locus in enumerate(loci, start=1):
        subject_genes = []
        for query_gene in locus:
            subject_gene = get_subject_gene(query_gene, data)
            if subject_gene:
                subject_genes.append(subject_gene)
        if subject_genes:
            subject_genes_str = ", ".join(subject_genes)
            print(f"Locus {i}: {subject_genes_str}")
            # Count core components in the current locus
            locus_core_count, locus_core_components = count_core_components_in_column(subject_genes)
            print(f"Arrangement of Locus {i}: {locus_core_count}")
        else:
            print(f"Locus {i}: No subject genes found")
if __name__ == "__main__":
    main()

