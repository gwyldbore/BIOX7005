import sys
import pandas as pd
import matplotlib.pyplot as plt

def main():
    input_files = snakemake.input

    # Initialize a DataFrame to store all replicates
    all_data = []

    # Load each replicate CSV and append to the DataFrame list
    for file in input_files:
        df = pd.read_csv(file)
        all_data.append(df)

    # Concatenate all replicate DataFrames into one
    combined_df = pd.concat(all_data)

    # sort by info column 
    combined_df = combined_df.sort_values('info')

    # find the points where the prediction changes
    prediction_changes = combined_df['overall_prediction'].shift() != combined_df['overall_prediction']

    # Filter rows where the prediction changes
    changed_df = combined_df[prediction_changes]

    # Create a frequency plot for 'num_mutation'
    plt.figure(figsize=(10, 6))
    changed_df['num_mutation'].value_counts().sort_index().plot(kind='bar', color='skyblue')

    # Add labels and title
    plt.title(f"Frequency of 'num_mutation' Changes")
    plt.xlabel("Number of Mutations")
    plt.ylabel("Frequency")

    # Save the plot
    plt.savefig(snakemake.output.graphs, bbox_inches='tight')

    


if __name__ == "__main__":
    main()
