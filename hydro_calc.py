#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 01:16:54 2024

@author: d18585wc
"""

import csv
import numpy as np
from Bio import SeqIO
from statistics import mean, stdev

def calculate_hydrophobicity(sequence, region_start, region_end):
    hydrophobic_aa = ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
    hydrophilic_aa = ['S', 'T', 'C', 'N', 'Q', 'D', 'E', 'K', 'R', 'H']
    amphiphatic_aa = ['G', 'P']
    
    region_sequence = sequence[region_start - 1:region_end]
    
    hydrophobic_count = sum(region_sequence.count(aa) for aa in hydrophobic_aa)
    hydrophilic_count = sum(region_sequence.count(aa) for aa in hydrophilic_aa)
    amphiphatic_count = sum(region_sequence.count(aa) for aa in amphiphatic_aa)
    
    total_aa_count = len(region_sequence)
    
    hydrophobic_percentage = (hydrophobic_count / total_aa_count) * 100
    hydrophilic_percentage = (hydrophilic_count / total_aa_count) * 100
    amphiphatic_percentage = (amphiphatic_count / total_aa_count) * 100
    
    return hydrophobic_percentage, hydrophilic_percentage, amphiphatic_percentage

def read_regions(region_file):
    regions_dict = {'N': {}, 'M': {}, 'C': {}}
    with open(region_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            sequence_id = row['sequence_header']
            try:
                regions_dict['N'][sequence_id] = (int(row['n_start']), int(row['n_end']))
                regions_dict['M'][sequence_id] = (int(row['m_start']), int(row['m_end']))
                regions_dict['C'][sequence_id] = (int(row['c_start']), int(row['c_end']))
            except ValueError as e:
                print(f"Warning: Skipping row with invalid data for sequence_header {sequence_id}: {e}")
                continue
    return regions_dict

def calculate_statistics(values):
    mean_value = mean(values)
    std_dev = stdev(values)
    sem = std_dev / np.sqrt(len(values))
    median = np.median(values)
    q1, q3 = np.percentile(values, [25, 75])
    return mean_value, std_dev, sem, median, q1, q3

def main():
    # File paths
    fasta_file = "all_seqs.txt"
    region_file = "hom_regions.csv"
    output_individual_file = "individual_region_stats.csv"
    output_overall_file = "overall_region_stats.csv"
    
    # Read region specifications
    regions_dict = read_regions(region_file)
    
    # Initialize storage for individual region statistics
    individual_stats = {'N': [], 'M': [], 'C': []}
    
    # Collect individual region stats
    with open(output_individual_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Sequence Header", "Region", "Hydrophobic (%)", "Hydrophilic (%)", "Amphiphatic (%)"])

        # Read sequences from FASTA file
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence_id = record.id
            sequence = str(record.seq)  # Convert sequence to string
            
            for region_type in ['N', 'M', 'C']:
                if sequence_id in regions_dict[region_type]:
                    region_start, region_end = regions_dict[region_type][sequence_id]
                    
                    if region_start < 1 or region_end > len(sequence):
                        print(f"Warning: Invalid region ({region_start}, {region_end}) for sequence {sequence_id}. Skipping region.")
                        continue
                    
                    hydrophobic_percentage, hydrophilic_percentage, amphiphatic_percentage = calculate_hydrophobicity(sequence, region_start, region_end)
                    
                    writer.writerow([sequence_id, region_type, hydrophobic_percentage, hydrophilic_percentage, amphiphatic_percentage])
                    
                    individual_stats[region_type].append({
                        'Sequence Header': sequence_id,
                        'Region': region_type,
                        'Hydrophobic (%)': hydrophobic_percentage,
                        'Hydrophilic (%)': hydrophilic_percentage,
                        'Amphiphatic (%)': amphiphatic_percentage
                    })
    
    # Calculate overall statistics for each region type (N, M, C)
    overall_stats = {}
    for region_type, stats in individual_stats.items():
        hydrophobic_values = [stat['Hydrophobic (%)'] for stat in stats]
        hydrophilic_values = [stat['Hydrophilic (%)'] for stat in stats]
        amphiphatic_values = [stat['Amphiphatic (%)'] for stat in stats]
        
        overall_stats[region_type] = {
            'Hydrophobic': calculate_statistics(hydrophobic_values),
            'Hydrophilic': calculate_statistics(hydrophilic_values),
            'Amphiphatic': calculate_statistics(amphiphatic_values)
        }
    
    # Save overall statistics to a CSV file
    with open(output_overall_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Region", "Statistic", "Mean", "Standard Deviation", "Standard Error of the Mean", "Median", "Q1", "Q3"])
        
        for region_type, region_stats in overall_stats.items():
            for aa_type, stats in region_stats.items():
                writer.writerow([region_type, f"{aa_type} (%)",
                                 stats[0],
                                 stats[1],
                                 stats[2],
                                 stats[3],
                                 stats[4],
                                 stats[5]])
    
if __name__ == "__main__":
    main()


   