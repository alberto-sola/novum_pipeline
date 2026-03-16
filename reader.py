# This script takes an RNAhybrid compact file and loads it into a pandas DataFrame.

from argparse import ArgumentParser  # Import the command-line argument parser.
import pandas as pd


COLUMNS = [
    "Gene",
    "Boh1",
    "miRNA",
    "miRNA_length",
    "Energy",
    "P_value",
    "Boh2",
    "miRNA_unmatches",
    "miRNA_matches",
    "target_matches",
    "target_unmatches",
]


def parse_args():
    parser = ArgumentParser(description="Load an RNAhybrid compact file into a pandas DataFrame.")  # Create the CLI parser.
    parser.add_argument("input_file", help="Path to the RNAhybrid compact TSV file.")  # Add the required input file argument.
    return parser.parse_args()  # Parse and return the arguments from the command line.


def load_predictions(input_file):
    df = pd.read_csv(input_file, sep=":", header=None)
    df.columns = COLUMNS  # Assign the predefined column names to the DataFrame.
    return df  # Return the populated DataFrame.


def main():
    args = parse_args()  # Read the file path passed by the user.
    df = load_predictions(args.input_file)  # Load the input file into a DataFrame.
    print(df.to_string(max_rows=30))  # Print the first five rows as a quick check.


if __name__ == "__main__":
    main()  # Run the script only when this file is executed directly.