
import random
from Bio.Seq import Seq
from tqdm import tqdm
from difflib import SequenceMatcher
from Bio import SeqIO

def jukes_cantor_mutate(sequence, mutation_rate, mutation_regions):
    nucleotides = ['A', 'C', 'G', 'T']
    mutated_sequence = list(sequence)  # Convert sequence to list for efficient mutation
    mutations = []
    for i, nucleotide in enumerate(sequence):
        if any(start <= i <= end for start, end in mutation_regions) and random.random() < mutation_rate:
            new_nucleotide = random.choice([nuc for nuc in nucleotides if nuc != nucleotide])
            mutated_sequence[i] = new_nucleotide
            mutations.append((i, nucleotide, new_nucleotide))
    return ''.join(mutated_sequence), mutations

def translate_sequence(nucleotide_sequence):
    seq = Seq(nucleotide_sequence)
    return seq.translate()

def check_stop_codons(amino_acid_sequence):
    return '*' in amino_acid_sequence

def calculate_charge(sequence, region_start, region_end, positive_threshold, negative_threshold):
    positive_aa = ['K', 'R', 'H']
    negative_aa = ['D', 'E']
    
    region_length = region_end - region_start + 1
    region_sequence = sequence[region_start - 1:region_end]
    
    positive_count = sum(region_sequence.count(aa) for aa in positive_aa)
    negative_count = sum(region_sequence.count(aa) for aa in negative_aa)
    
    total_aa_count = len(region_sequence)
    
    positive_percentage = (positive_count / total_aa_count) * 100
    negative_percentage = (negative_count / total_aa_count) * 100
    
    return positive_percentage, negative_percentage


# Read nucleotide sequence from file
def read_sequence_from_file(filename):
    with open(filename, 'r') as file:
        sequence = file.read().strip().lower()  # Read sequence in lowercase
    return sequence

# Input filename containing nucleotide sequence
filename = "futsch_dna.txt"

# Output filename for mutation logs
output_filename = "Oops_Logs_and_Similarity.txt"


# Read nucleotide sequence from file
nucleotide_sequence = read_sequence_from_file(filename)

# Mutation rate for Jukes-Cantor model (per nucleotide per generation)
mutation_rate = 4.5e-5

# Drosophila estimated mutation rate = 2.8e-9

# Number of generations
num_generations = 200000

# Initialize variables
current_sequence = nucleotide_sequence.upper()  # Start with uppercase for mutation compatibility
total_mutations = 0
stop_codon_logs = []
mutation_logs = []
charge_failures = []


# Define the two predefined regions for mutation
mutation_regions = [(1782, 16146)]  # Example: Mutate regions from 100 to 199 and 500 to 599

# Define the three predefined regions for charge
charge_regions = [(741, 929), (930, 5340), (5341, 5410)]

# Define positive and negative charge thresholds for each charge region
charge_thresholds = {
    1: {'positive': (22.87, 26.47), 'negative': (4.88, 8.88)},
    2: {'positive': (13.76, 17.76), 'negative': (19.51, 23.51)},
    3: {'positive': (18.00, 22.00), 'negative': (9.43, 13.43)},
}

# Move the translation function outside the mutation loop
original_amino_acid_sequence = translate_sequence(nucleotide_sequence)

# Move the file reading outside the loop
nucleotide_sequence = read_sequence_from_file(filename)


# Progress bar
with tqdm(total=num_generations, desc="Mutating sequence", dynamic_ncols=True) as pbar:
    for generation in range(num_generations):
        # Store total mutations before mutation
        total_mutations_before_mutation = total_mutations
        
        # Store the last valid mutated sequence before any reversion takes place
        last_valid_sequence = current_sequence
        
        # Calculate the number of mutations for this generation based on the mutation rate
        num_mutations = int(len(current_sequence) * mutation_rate)
        
        # Perform mutation using previous sequence
        mutated_sequence, mutations = jukes_cantor_mutate(current_sequence, mutation_rate, mutation_regions)
        
        # Calculate total mutations after mutation
        total_mutations += len(mutations)
        
        # Perform further operations only if mutations were introduced
        if total_mutations > total_mutations_before_mutation:
            # Translate the mutated sequence into amino acids
            mutated_amino_acid_sequence = translate_sequence(mutated_sequence)
            
            # Check for stop codons in the translated amino acid sequence
            stop_codons_present = check_stop_codons(mutated_amino_acid_sequence)
            
            # Perform charge threshold check for each charge region
            charge_threshold_passed = True
            for i, (region_start, region_end) in enumerate(charge_regions, start=1):
                positive_threshold, negative_threshold = charge_thresholds[i]['positive'], charge_thresholds[i]['negative']
                positive_percentage, negative_percentage = calculate_charge(mutated_amino_acid_sequence, region_start, region_end, positive_threshold, negative_threshold)
                if not (positive_threshold[0] <= positive_percentage <= positive_threshold[1] and negative_threshold[0] <= negative_percentage <= negative_threshold[1]):
                    charge_threshold_passed = False
                    charge_failures.append((generation, i, positive_percentage, negative_percentage))
            
          

            
            # Revert to the previous mutated sequence if stop codons are present, charge thresholds are not met, or hydrophobicity thresholds are not met
            if stop_codons_present or not charge_threshold_passed:
                current_sequence = last_valid_sequence  # Revert to the last valid mutated sequence
            else:
                current_sequence = mutated_sequence.upper()  # Store mutated sequence in uppercase
            
            # Record mutation logs if mutations occurred
            if mutations:
                mutation_logs.append((generation, mutations))
        
        pbar.update(1)  # Update progress bar
        
# Calculate final mutated amino acid sequence
final_mutated_amino_acid_sequence = translate_sequence(current_sequence)

#Define the original amino acid sequence
first_protein_sequence = translate_sequence(nucleotide_sequence)

# Prepare charge calculations for each region based on the final mutated amino acid sequence
charge_results = {}
for i, (region_start, region_end) in enumerate(charge_regions, start=1):
    positive_threshold, negative_threshold = charge_thresholds[i]['positive'], charge_thresholds[i]['negative']
    positive_percentage, negative_percentage = calculate_charge(final_mutated_amino_acid_sequence, region_start, region_end, positive_threshold, negative_threshold)
    charge_results[i] = {'positive_percentage': positive_percentage, 'negative_percentage': negative_percentage}

# Output results to file
with open(output_filename, 'w') as outfile:
    # Write mutation logs
    outfile.write("\nMutation Logs:\n")
    for generation, mutations in mutation_logs:
        if mutations:  # Only write if there were mutations in this generation
            outfile.write(f"Generation {generation + 1}:\n")
            for mutation in mutations:
                outfile.write(f"Position: {mutation[0]}, Original: {mutation[1]}, Mutated: {mutation[2]}\n")

    # Write charge failure logs
    outfile.write("\nCharge Failure Logs:\n")
    for generation, region, positive_percentage, negative_percentage in charge_failures:
        outfile.write(f"Generation {generation + 1}: Charge threshold not met for region {region}, Positive: {positive_percentage}%, Negative: {negative_percentage}%\n")

    # Write charge calculations for each region based on the final mutated amino acid sequence
    outfile.write("Charge Calculations for Each Region:\n")
    for region, charge_data in charge_results.items():
        outfile.write(f"Region {region}:\n")
        outfile.write(f"Positive Percentage: {charge_data['positive_percentage']:.2f}%\n")
        outfile.write(f"Negative Percentage: {charge_data['negative_percentage']:.2f}%\n")
        outfile.write("\n")

# Output results
print("\nOriginal Nucleotide Sequence:\n", nucleotide_sequence)
print("\nOriginal Amino Acid Sequence:\n", first_protein_sequence)
print("\nFinal Mutated Nucleotide Sequence:\n", current_sequence)
print("\nMutated Amino Acid Sequence:\n", final_mutated_amino_acid_sequence)
print("\nTotal Mutations:\n", total_mutations)
