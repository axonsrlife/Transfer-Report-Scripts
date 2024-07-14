import os
import subprocess
import urllib.request

# Function to check if a file is a valid FASTA file
def is_valid_fasta(file_path):
    with open(file_path, 'rb') as f:
        first_line = f.readline().strip()
        return first_line.startswith(b'>')

# Step 1: Check if the database exists, if not download and make it
def prepare_database(genome_url, genome_file, database_name):
    if not os.path.exists(genome_file):
        print("Downloading Drosophila melanogaster genome...")
        urllib.request.urlretrieve(genome_url, genome_file)
        print("Genome downloaded.")

    if not os.path.exists(database_name + ".phr"):
        print("Creating BLAST database...")
        # Check if the input genome file is a valid FASTA file
        if is_valid_fasta(genome_file):
            subprocess.run(['makeblastdb', '-in', genome_file, '-dbtype', 'prot', '-out', database_name])
            print("BLAST database created.")
        else:
            print("Error: Input file is not a valid FASTA file.")

# Step 2: Run BLASTp
def run_blastp(query_sequence, database_name, output_file):
    print("Running BLASTp...")
    if os.path.exists(database_name + ".phr"):
        subprocess.run(['blastp', '-query', query_sequence, '-db', database_name, '-out', output_file, '-evalue', '1e-5', '-outfmt', '6'])
        print("BLASTp finished. Results saved to:", output_file)
    else:
        print("Error: BLAST database not found. Please check if the database was created successfully.")

# Step 3: Extract homologous regions to futsch protein
def extract_homologous_regions(blast_result_file, futsch_protein_id, output_dir):
    print("Extracting homologous regions to futsch protein...")
    homologous_regions = []
    with open(blast_result_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            query_id = fields[0]
            subject_id = fields[1]
            if subject_id == futsch_protein_id:
                homologous_regions.append((query_id, int(fields[6]), int(fields[7])))
    print("Homologous regions extracted.")

    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Write homologous regions to file
    with open(os.path.join(output_dir, "homologous_regions.txt"), 'w') as out_f:
        for region in homologous_regions:
            out_f.write(f"{region[0]}\t{region[1]}\t{region[2]}\n")

# Main function
def main():
    # Define paths and parameters
    genome_url = 'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001215.4/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT'
    genome_file = 'drosophila_melanogaster_genome.fasta'
    database_name = 'drosophila_melanogaster_db'
    query_sequence_file = 'seqs_test2.txt'  # Change this to the path of your query sequence file
    output_file = 'blastp_results.txt'
    futsch_protein_id = 'Q9W596'
    output_dir = 'homologous_regions'

    # Step 1: Prepare the database
    prepare_database(genome_url, genome_file, database_name)

    # Step 2: Run BLASTp
    run_blastp(query_sequence_file, database_name, output_file)

    # Step 3: Extract homologous regions
    extract_homologous_regions(output_file, futsch_protein_id, output_dir)

if __name__ == "__main__":
    main()
