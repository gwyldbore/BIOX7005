import sys
import pandas as pd
import matplotlib.pyplot as plt

def main():
    input_files = snakemake.input

    all_mutation_counts = []

    for file in input_files:
        # Load the data from the file
        df = pd.read_csv(file)

        # Sort by 'info' to ensure the transitions are in correct order
        df = df.sort_values('info')

        # Identify the starting value (first prediction in this file)
        starting_value = df['overall_prediction'].iloc[0]

        # Define the other categories to track (excluding the starting value)
        categories_to_track = ["NR1", "NR1-like", "NR4-like", "NR4"]
        categories_to_track.remove(starting_value)

        # Track the first transition to each non-starting category
        for category in categories_to_track:
            # Identify where the prediction changes to the target category
            transition = df[(df['overall_prediction'].shift() != df['overall_prediction']) & 
                            (df['overall_prediction'] == category)]

            # If there is at least one such transition, store the first one
            if not transition.empty:
                first_transition = transition.iloc[0]
                all_mutation_counts.append(first_transition['num_mutation'])

    # Create a frequency plot from the aggregated mutation counts
    plt.figure(figsize=(10, 6))
    pd.Series(all_mutation_counts).value_counts().sort_index().plot(kind='bar', color='skyblue')

    # Add labels and title
    plt.title(f"Frequency of 'num_mutation' Changes")
    plt.xlabel("Number of Mutations")
    plt.ylabel("Frequency")

    # Save the plot
    plt.savefig(snakemake.output.graphs, bbox_inches='tight')

    


if __name__ == "__main__":
    main()
