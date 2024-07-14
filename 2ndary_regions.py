import csv

def read_fasta(filename):
    sequences = {}
    with open(filename, "r") as file:
        sequence_id = None
        sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if sequence_id:
                    sequences[sequence_id] = sequence
                sequence_id = line[1:]
                sequence = ""
            else:
                if sequence_id:
                    sequence += line
    if sequence_id:
        sequences[sequence_id] = sequence
    return sequences

def extract_structure(sequence):
    half_length = len(sequence) // 2
    amino_acids = sequence[:half_length]
    structure = sequence[half_length:]
    return amino_acids, structure

def generate_annotation(amino_acids, structure):
    annotation = ""
    for aa, ss in zip(amino_acids, structure):
        if ss == "H":
            annotation += "H"
        elif ss == "E":
            annotation += "E"
        elif ss == "C":
            annotation += "-"
        elif ss == "T":
            annotation += "T"
    return annotation

def find_structure_regions(annotation):
    regions = []
    current_region_start = None
    current_region_type = None
    for i, ss in enumerate(annotation):
        if ss != '-':
            if current_region_start is None:
                current_region_start = i
                current_region_type = ss
            elif ss != current_region_type:
                current_region_end = i - 1
                if current_region_end - current_region_start >= 10:
                    regions.append((current_region_start, current_region_end, current_region_type))
                current_region_start = i
                current_region_type = ss
    if current_region_start is not None:
        current_region_end = len(annotation) - 1
        if current_region_end - current_region_start >= 10:
            regions.append((current_region_start, current_region_end, current_region_type))
    return regions

if __name__ == "__main__":
    filename = "./s4pred/output3.fas"  # Change this to your output FASTA file
    sequences = read_fasta(filename)
    with open("Structure_Regions.csv", "w", newline='') as output_file:
        csv_writer = csv.writer(output_file)
        csv_writer.writerow(["Sequence ID", "Start", "End", "Type"])
        
        for sequence_id, sequence in sequences.items():
            amino_acids, structure = extract_structure(sequence)
            annotation = generate_annotation(amino_acids, structure)
            regions = find_structure_regions(annotation)
            
            for start, end, region_type in regions:
                csv_writer.writerow([sequence_id, start, end, region_type])
