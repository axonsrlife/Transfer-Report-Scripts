import matplotlib.pyplot as plt
import pandas as pd

# Function to read domain data from CSV
def read_domain_data_from_csv(file_path):
    df = pd.read_csv(file_path)
    domains = {}
    for _, row in df.iterrows():
        query_id = row['Query_ID']
        start = row['Subject_Start']
        end = row['Subject_End']
        if query_id not in domains:
            domains[query_id] = []
        domains[query_id].append((start, end))
    return domains

# Function to plot domain diagram
def plot_domain_diagram(domains):
    fig, ax = plt.subplots()
    for i, query_id in enumerate(domains):
        for start, end in domains[query_id]:
            # Plot each domain as a rectangle on the same line for the same Query ID
            ax.add_patch(plt.Rectangle((start, i), end - start, 0.8, color='blue', alpha=0.5))
    ax.set_xlabel('Amino Acid Position')
    ax.set_ylabel('Query ID')
    ax.set_title('Homologous Domain Diagram')
    
    # Adjust y-axis limits and ticks
    ax.set_ylim(-1, len(domains))
    ax.set_yticks(range(len(domains)))
    ax.set_yticklabels(list(domains.keys()))
    
    # Adjust x-axis limits
    x_max = max(max(end for _, end in domains[query_id]) for query_id in domains)
    ax.set_xlim(0, x_max)
    
    plt.grid(True)  # Add gridlines for better visualization
    plt.show()

# Path to the CSV file containing domain data
csv_file_path = 'futsch_homology.csv'

# Read domain data from the CSV file
domains = read_domain_data_from_csv(csv_file_path)

# Plot the domain diagram
plot_domain_diagram(domains)
