import sys
import pandas as pd
import matplotlib.pyplot as plt

def main():
    input_files = snakemake.input

    # Initialize a dictionary to store mutation counts for each target category
    mutation_counts = {
        "NR1": [],
        "NR1-like": [],
        "NR4-like": [],
        "NR4": []
    }

    # Loop through each input file and extract transitions
    for file in input_files:
        # Load the data from the file
        df = pd.read_csv(file)

        # Identify the starting value (first prediction in this file)
        starting_value = df['overall_prediction'].iloc[0]

        # Define the other categories to track (excluding the starting value)
        categories_to_track = ["NR1", "NR1-like", "NR4-like", "NR4"]
        categories_to_track.remove(starting_value)

        if starting_value == 'NR4':
            categories_to_track.reverse()

        # Track the first transition to each non-starting category
        for category in categories_to_track:
            # Identify where the prediction changes to the target category
            transition = df[(df['overall_prediction'].shift() != df['overall_prediction']) & 
                            (df['overall_prediction'] == category)]

            # If there is at least one such transition, store the first one
            if not transition.empty:
                first_transition = transition.iloc[0]
                mutation_counts[category].append(first_transition['num_mutation'])

    # Create a grid of 3 subplots (one for each non-starting category)
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

    # print(mutation_counts.items())
    # Plot the frequency for each target category in a separate subplot
    for ax, (category, counts) in zip(axes, categories_to_track):
        print(f'category is: {category}')
        if counts:  # Only plot if there are relevant transitions
            pd.Series(counts).value_counts().sort_index().plot(kind='bar', color='skyblue', ax=ax)
        ax.set_title(f"Transitions to {category}")
        ax.set_xlabel("Number of Mutations")
        ax.set_ylabel("Frequency")

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save the plot
    plt.savefig(snakemake.output.graphs, bbox_inches='tight')

    


if __name__ == "__main__":
    main()
