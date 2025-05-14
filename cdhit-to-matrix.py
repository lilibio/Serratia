import re
import sys

def parse_clusters(input_file):
    clusters = {}
    current_cluster = None

    with open(input_file, 'r') as file:
        for line in file:
            if line.startswith('>Cluster'):
                current_cluster = line.strip()
                clusters[current_cluster] = []
            elif line.strip():
                gene_match = re.search(r'>(\S+)', line)
                if gene_match:
                    gene_name = gene_match.group(1).rstrip("...")
                    clusters[current_cluster].append(gene_name)

    return clusters

def format_clusters(clusters):
    formatted_clusters = []
    for cluster, genes in clusters.items():
        cluster_name = f"Cluster_{cluster.split()[-1]}"
        formatted_cluster = [cluster_name] + genes
        formatted_clusters.append('\t'.join(formatted_cluster))

    return formatted_clusters

def write_output(output_file, formatted_clusters):
    with open(output_file, 'w') as file:
        file.write('\n'.join(formatted_clusters))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file output_file")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    clusters = parse_clusters(input_file)
    formatted_clusters = format_clusters(clusters)
    write_output(output_file, formatted_clusters)

    print(f"Formatted clusters written to {output_file}")
