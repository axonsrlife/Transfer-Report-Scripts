from Bio import SeqIO
from Bio.Seq import Seq

def calculate_charge(sequence, region_start, region_end):
    positive_aa = ['K', 'R', 'H']
    negative_aa = ['D', 'E']
    
    region_sequence = sequence[region_start - 1:region_end]
    
    positive_count = sum(region_sequence.count(aa) for aa in positive_aa)
    negative_count = sum(region_sequence.count(aa) for aa in negative_aa)
    
    total_aa_count = len(region_sequence)
    
    if total_aa_count == 0:
        return 0, 0
    
    positive_percentage = (positive_count / total_aa_count) * 100
    negative_percentage = (negative_count / total_aa_count) * 100
    
    return positive_percentage, negative_percentage

def main():
    # Read the FASTA file containing nucleotide sequence
    fasta_file = "futsch_mdna.fasta"
    record = SeqIO.read(fasta_file, "fasta")
    
    # Translate the nucleotide sequence to protein sequence
    nucleotide_sequence = record.seq
    protein_sequence = nucleotide_sequence.translate()
    
    # Define the five predefined regions
    regions = [(1, 740), (741, 929), (930, 5340), (5341, 5410), (5411, 5495)]
    
    print("Region\tPositive Charge (%)\tNegative Charge (%)")
    for i, region in enumerate(regions, start=1):
        region_start, region_end = region
        
        # Calculate charge percentages for the region
        positive_percentage, negative_percentage = calculate_charge(protein_sequence, region_start, region_end)
        
        print(f"Region {i}\t{positive_percentage:.2f}\t\t\t{negative_percentage:.2f}")

if __name__ == "__main__":
    main()
