import csv
from statistics import mean, stdev
from math import sqrt

def read_fasta(filename):
    sequences = {}
    with open(filename, "r") as file:
        sequence_id = None
        sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if sequence_id:
                    sequences[sequence_id.lower()] = sequence  # Convert to lowercase
                sequence_id = line[1:]
                sequence = ""
            else:
                if sequence_id:
                    sequence += line
    if sequence_id:
        sequences[sequence_id.lower()] = sequence  # Convert to lowercase
    return sequences

def extract_structure(sequence):
    half_length = len(sequence) // 2
    amino_acids = sequence[:half_length]
    structure = sequence[half_length:]
    return amino_acids, structure

def read_defined_regions_csv(filename):
    defined_regions = {}
    with open(filename, "r") as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            sequence_id = row['sequence_header'].lower()  # Convert to lowercase
            # Adjust for zero-based indexing
            n_start = int(row['n_start']) - 1 if row['n_start'] else None
            n_end = int(row['n_end']) - 1 if row['n_end'] else None
            m_start = int(row['m_start']) - 1 if row['m_start'] else None
            m_end = int(row['m_end']) - 1 if row['m_end'] else None
            c_start = int(row['c_start']) - 1 if row['c_start'] else None
            c_end = int(row['c_end']) - 1 if row['c_end'] else None

            if sequence_id not in defined_regions:
                defined_regions[sequence_id] = []

            if n_start is not None and n_end is not None:
                defined_regions[sequence_id].append({
                    'region_name': 'N-terminal',
                    'start': n_start,
                    'end': n_end
                })
            if m_start is not None and m_end is not None:
                defined_regions[sequence_id].append({
                    'region_name': 'M-region',
                    'start': m_start,
                    'end': m_end
                })
            if c_start is not None and c_end is not None:
                defined_regions[sequence_id].append({
                    'region_name': 'C-terminal',
                    'start': c_start,
                    'end': c_end
                })
    return defined_regions

def calculate_structure_percentage(amino_acids, structure, defined_regions):
    total_secondary_length = sum(1 for s in structure if s in {'H', 'E'})
    percentages = []

    max_secondary_percentage = 0
    max_total_secondary_percentage = 0

    for region in defined_regions:
        start = region['start']
        end = region['end']
        region_name = region['region_name']
        region_length = end - start + 1
        if region_length <= 0:
            continue

        secondary_count = 0
        for i in range(start, end + 1):
            if i >= len(structure):  # Check if the index is within the structure length
                print(f"Index {i} out of range for structure of length {len(structure)}")
                continue
            if structure[i] in {'H', 'E'}:
                secondary_count += 1

        percentage_in_region = (secondary_count / region_length) * 100
        percentage_of_total_secondary = (secondary_count / total_secondary_length) * 100 if total_secondary_length > 0 else 0

        # Track max percentages for normalization
        max_secondary_percentage = max(max_secondary_percentage, percentage_in_region)
        max_total_secondary_percentage = max(max_total_secondary_percentage, percentage_of_total_secondary)

        percentages.append({
            "Region": region_name,
            "Start": start + 1,  # Convert back to one-based for output
            "End": end + 1,      # Convert back to one-based for output
            "Secondary_Percentage": percentage_in_region,
            "Percentage_of_Total_Secondary": percentage_of_total_secondary
        })

    for region_info in percentages:
        percentage_in_region = region_info['Secondary_Percentage']
        percentage_of_total_secondary = region_info['Percentage_of_Total_Secondary']

        # Combined Structural Score
        combined_score = (percentage_in_region + percentage_of_total_secondary) / 2

        # Normalized Score
        normalized_secondary = (percentage_in_region / max_secondary_percentage) if max_secondary_percentage else 0
        normalized_total_secondary = (percentage_of_total_secondary / max_total_secondary_percentage) if max_total_secondary_percentage else 0
        normalized_score = normalized_secondary + normalized_total_secondary

        # Ratio
        ratio = (percentage_in_region / percentage_of_total_secondary) if percentage_of_total_secondary else 0

        region_info.update({
            "Combined_Score": combined_score,
            "Normalized_Score": normalized_score,
            "Ratio": ratio
        })

    return percentages

def calculate_descriptive_statistics(percentages):
    stats = {}
    region_types = ['N-terminal', 'M-region', 'C-terminal']
    metrics = ["Secondary_Percentage", "Percentage_of_Total_Secondary", "Combined_Score", "Normalized_Score", "Ratio"]

    for region in region_types:
        stats[region] = {metric: [] for metric in metrics}

    for entry in percentages:
        for metric in metrics:
            stats[entry["Region"]][metric].append(entry[metric])

    descriptive_stats = []

    for region in region_types:
        for metric in metrics:
            if stats[region][metric]:
                mean_value = mean(stats[region][metric])
                std_dev = stdev(stats[region][metric]) if len(stats[region][metric]) > 1 else 0
                std_err = std_dev / sqrt(len(stats[region][metric])) if len(stats[region][metric]) > 1 else 0
                min_value = min(stats[region][metric])
                max_value = max(stats[region][metric])
                descriptive_stats.append({
                    "Region": region,
                    "Metric": metric,
                    "Mean": mean_value,
                    "Std_Dev": std_dev,
                    "Std_Err": std_err,
                    "Min": min_value,
                    "Max": max_value
                })

    return descriptive_stats

def print_sequence_lengths(sequences):
    for seq_id, sequence in sequences.items():
        amino_acids, structure = extract_structure(sequence)
        print(f"Sequence ID: {seq_id}, Length: {len(sequence)}, Amino Acids Length: {len(amino_acids)}, Structure Length: {len(structure)}")

if __name__ == "__main__":
    fasta_filename = "./s4pred/output3.fas"  # Path to the s4pred output FASTA file
    defined_regions_csv_filename = "hom_regions.csv"  # Path to the defined regions CSV
    output_filename = "Secondary_Structure_Percentages.csv"  # Output CSV file
    stats_output_filename = "Secondary_Structure_Statistics.csv"  # Statistics CSV file

    sequences = read_fasta(fasta_filename)
    defined_regions = read_defined_regions_csv(defined_regions_csv_filename)

    # Print lengths of sequences for debugging
    print_sequence_lengths(sequences)

    all_percentages = []

    with open(output_filename, "w", newline='') as output_file:
        csv_writer = csv.writer(output_file)
        csv_writer.writerow([
            "Sequence ID",
            "Region",
            "Region Start",
            "Region End",
            "Secondary Structure Percentage",
            "Percentage of Total Secondary Structure",
            "Combined Structural Score",
            "Normalized Score",
            "Ratio"
        ])

        for sequence_id, sequence in sequences.items():
            if sequence_id not in defined_regions:
                print(f"No defined regions for sequence {sequence_id}")
                continue

            amino_acids, structure = extract_structure(sequence)
            percentages = calculate_structure_percentage(amino_acids, structure, defined_regions[sequence_id])
            all_percentages.extend(percentages)

            for region_info in percentages:
                csv_writer.writerow([
                    sequence_id,
                    region_info["Region"],
                    region_info["Start"],
                    region_info["End"],
                    f"{region_info['Secondary_Percentage']:.2f}",
                    f"{region_info['Percentage_of_Total_Secondary']:.2f}",
                    f"{region_info['Combined_Score']:.2f}",
                    f"{region_info['Normalized_Score']:.2f}",
                    f"{region_info['Ratio']:.2f}"
                ])

    descriptive_stats = calculate_descriptive_statistics(all_percentages)

    with open(stats_output_filename, "w", newline='') as stats_file:
        csv_writer = csv.writer(stats_file)
        csv_writer.writerow(["Region", "Metric", "Mean", "Std_Dev", "Std_Err", "Min", "Max"])

        for stat in descriptive_stats:
            csv_writer.writerow([
                stat["Region"],
                stat["Metric"],
                f"{stat['Mean']:.2f}",
                f"{stat['Std_Dev']:.2f}",
                f"{stat['Std_Err']:.2f}",
                f"{stat['Min']:.2f}",
                f"{stat['Max']:.2f}"
            ])
