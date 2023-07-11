from Bio import SeqIO
import csv
import pandas as pd
import matplotlib.pyplot as plt


def CleanFile():
    fasta_file = "data/total_msa.fasta"
    csv_file = "data/clean.csv"

    # Open the CSV file for writing
    with open(csv_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Taxon Name','Length','City', 'Sequence'])  # Write the header

        # Iterate over each sequence in the FASTA file
        for record in SeqIO.parse(fasta_file, "fasta"):
            taxon_name = record.id
            sequence = str(record.seq)

            # Write the taxon name and sequence to the CSV file
            writer.writerow([taxon_name,len(sequence),taxon_name.split("-")[0] ,sequence])

    print("CSV file created successfully.")
def MakePlot():
    csv_file = "data/clean.csv"

    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(csv_file)

    # Extract the city names from the taxon names
    df['City'] = df['Taxon Name'].str.split('-', expand=True)[0]

    # Count the frequency of each city
    city_counts = df['City'].value_counts()

    # Plot the frequency of city names
    plt.figure(figsize=(10, 6))
    city_counts.plot(kind='bar')
    plt.xlabel('City')
    plt.ylabel('Frequency')
    plt.title('Frequency of City Names')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()
def findProfile():
    from Bio import SeqIO
    import pandas as pd

    fasta_file = "data/total_msa.fasta"

    # Initialize a dictionary to store the frequency of mutants at each position
    position_counts = {}

    # Iterate over each sequence in the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)

        # Iterate over each position in the sequence
        for i, nucleotide in enumerate(sequence):
            # Skip any gaps or missing data
            if nucleotide == "-" or nucleotide == "?":
                continue

            # Update the frequency count for the current position and nucleotide
            if i in position_counts:
                if nucleotide in position_counts[i]:
                    position_counts[i][nucleotide] += 1
                else:
                    position_counts[i][nucleotide] = 1
            else:
                position_counts[i] = {nucleotide: 1}

    # Convert the frequency counts into a pandas DataFrame
    df = pd.DataFrame(position_counts)

    # Fill any missing positions with zero frequency
    df.fillna(0, inplace=True)
    df['Total'] = df.sum(axis=1)
    df.drop('Total', axis=1, inplace=True)
    df.sort_index(axis=1, inplace=True)

    df.to_csv('profile.csv', index_label='Position')

    print("Frequency of mutants at each position saved successfully.")

def MottifMutant():
    from Bio import SeqIO
    import pandas as pd

    fasta_file = "data/total_msa.fasta"
    window_size = 10

    # Initialize lists to store window positions and mutant frequencies
    positions = []
    mutant_frequencies = []

    # Iterate over each sequence in the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)

        # Iterate over each window in the sequence
        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i + window_size]

            # Calculate the frequency of mutants in the window
            mutant_count = sum(1 for nucleotide in window if nucleotide.islower())
            mutant_frequency = mutant_count / window_size

            positions.append(i)
            mutant_frequencies.append(mutant_frequency)

    # Convert the data into a pandas DataFrame
    df = pd.DataFrame({'Position': positions, 'Mutant Frequency': mutant_frequencies})
    max_region = df.loc[df['Mutant Frequency'].idxmax()]

    print("Region with the highest mutant frequency:")
    print(max_region)

import numpy as np

def find_maximum_mutant_region(fasta_file):
    sequences = []
    df = pd.read_csv(fasta_file)
    data = df.to_numpy()
    data = data[:,3]
    for row in data:
        sequences.append(list(row))
    alignment = np.array(sequences)
    alignment_length = alignment.shape[1]

    mutant_region = np.zeros(alignment_length)
    mutant_region_dict = []
    for i in range(len(mutant_region)):
        mutant_region_dict.append(dict())

    wuhan = alignment[0,:]
    for i in range(1,alignment.shape[0]):
        for j in range(alignment.shape[1]):
            if (alignment[i,j]!=wuhan[j]):
                mutant_region[j]+=1
                nucleutid = alignment[i,j]
                try:
                    mutant_region_dict[j][nucleutid]+=1
                except Exception as e:
                    mutant_region_dict[j][nucleutid] = 1

    for i, item in enumerate(mutant_region_dict):
        if len(item)>=3:
            print(i)
    x = range(1,len(wuhan)+1)
    plt.figure(figsize=(10, 6))
    plt.bar(x, mutant_region)
    plt.xlabel('Position')
    plt.ylabel('Frequency')
    plt.title('Frequency Mutant')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()


# Example usage
fasta_file = 'data/clean.csv'
find_maximum_mutant_region(fasta_file)



#CleanFile()
#MakePlot()
#findProfile()
#MottifMutant()