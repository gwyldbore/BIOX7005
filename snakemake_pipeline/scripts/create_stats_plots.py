import sys
import pandas as pd
import matplotlib.pyplot as plt
from ast import literal_eval

def extract_mutation_counts(input_files) -> dict:
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
            mutation_counts.pop('NR4', '')
        elif starting_value == 'NR1':
            mutation_counts.pop('NR1', '')

        # Track the first transition to each non-starting category
        for category in categories_to_track:
            # Identify where the prediction changes to the target category
            transition = df[(df['overall_prediction'].shift() != df['overall_prediction']) & (df['overall_prediction'] == category)]

            # If there is at least one such transition, store the first one
            if not transition.empty:
                first_transition = transition.iloc[0]
                mutation_counts[category].append(first_transition['num_mutation'])

    # Reorder mutation counts based on the valid categories order
    ordered_counts = {cat: mutation_counts[cat] for cat in categories_to_track if cat in mutation_counts}

    return ordered_counts


def extract_mutated_positions(input_files):
    all_categories = ["NR1", "NR1-like", "NR4-like", "NR4"]
    mutated_positions = {category: [] for category in all_categories}
    max_sequence_length = 0

        # Loop through each input file and extract transitions
    for file in input_files:
        # Load the data from the file
        df = pd.read_csv(file)

        df['mutated_positions'] = df['mutated_positions'].apply(lambda x: literal_eval(str(x)))

        df['sequence_length'] = df['sequence'].apply(lambda x: len(str(x)) if pd.notna(x) else 0)
        max_sequence_length = max(max_sequence_length, df['sequence_length'].max())

        # Identify the starting value (first prediction in this file)
        starting_value = df['overall_prediction'].iloc[0]

        # Define the other categories to track (excluding the starting value)
        categories_to_track = ["NR1", "NR1-like", "NR4-like", "NR4"]
        categories_to_track.remove(starting_value)
        
        if starting_value == 'NR4':
            categories_to_track.reverse()
            mutated_positions.pop('NR4', '')
        elif starting_value == 'NR1':
            mutated_positions.pop('NR1', '')

        for category in categories_to_track:
            transition = df[(df['overall_prediction'].shift() != df['overall_prediction']) &
                            (df['overall_prediction'] == category)]
            if not transition.empty:
                first_transition = transition.iloc[0]
                mutated_positions[category].extend(first_transition['mutated_positions'])

        ordered_positions = {cat: mutated_positions[cat] for cat in categories_to_track if mutated_positions[cat]}
        return ordered_positions, max_sequence_length



def plot_num_mutations(ordered_counts, output_path):
    # Create a grid of 3 subplots (one for each non-starting category)
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

    # Plot the frequency for each target category in a separate subplot
    for ax, (category, counts) in zip(axes, ordered_counts.items()):
        if counts:  # Only plot if there are relevant transitions
            pd.Series(counts).value_counts().sort_index().plot(kind='bar', color='skyblue', ax=ax)
        ax.set_title(f"Transitions to {category}")
        ax.set_xlabel("Number of Mutations")
        ax.set_ylabel("Frequency")

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_path, bbox_inches='tight')


def plot_mutated_positions(ordered_positions, sequence_length, output_path):
    # fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

    # for ax, (category, positions) in zip(axes, ordered_positions.items()):
    #     # Calculate frequency of mutated positions
    #     counts = pd.Series(positions).value_counts().sort_index()

    #     # Plot the frequency of mutated positions as bars
    #     counts.plot(kind='bar', color='skyblue', ax=ax)

    #     # Only plot every 10th tick if the sequence is too long
    #     if sequence_length > 50:
    #         tick_spacing = 10
    #     else:
    #         tick_spacing = 1

    #     ax.set_xticks(range(0, sequence_length + 1, tick_spacing))
    #     ax.set_xticklabels(range(0, sequence_length + 1, tick_spacing), rotation=90)

    #     # Set plot title and labels
    #     ax.set_title(f"Mutated Positions to {category}")
    #     ax.set_xlabel("Sequence Position")
    #     ax.set_ylabel("Frequency")

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

    print(f'mutated positions list being plotted: {[key, value.sorted() for key, value in ordered_positions.items()]}')

    for ax, (category, positions) in zip(axes, ordered_positions.items()):
        # Calculate frequency of each mutation position
        counts = pd.Series(positions).value_counts().sort_index()

        # Plot the mutation positions as bars with the correct x-axis values
        ax.bar(counts.index, counts.values, color='skyblue')

        # Set x-axis range to cover the full sequence length
        ax.set_xlim(0, sequence_length)

        # Set x-ticks for every position in the range of the sequence
        ax.set_xticks(range(0, sequence_length + 1, max(1, sequence_length // 20)))  # Adjust tick density

        # Rotate x-tick labels for readability
        ax.set_xticklabels(range(0, sequence_length + 1, max(1, sequence_length // 20)), rotation=90)

        # Set plot title and labels
        ax.set_title(f"Mutated Positions to {category}")
        ax.set_xlabel("Sequence Position")
        ax.set_ylabel("Frequency")

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_path, bbox_inches='tight')




def main():
    input_files = snakemake.input

    ordered_counts = extract_mutation_counts(input_files)
    plot_num_mutations(ordered_counts, snakemake.output.mutation_graphs)

    ordered_positions, sequence_length = extract_mutated_positions(input_files)
    plot_mutated_positions(ordered_positions, sequence_length, snakemake.output.position_graphs)

    


                
    

    


if __name__ == "__main__":
    main()
