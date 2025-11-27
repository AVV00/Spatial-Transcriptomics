import argparse
import pandas as pd
import glob

def combine_files(file_pattern, output_file):
    # Get a list of file paths matching the pattern
    file_paths = glob.glob(file_pattern)

    # Initialize an empty DataFrame to store the combined data
    combined_data = pd.DataFrame()

    # Iterate over the file paths and read each TSV file
    for file_path in file_paths:
        # Read the TSV file into a DataFrame
        data = pd.read_csv(file_path, delimiter='\t', index_col=0, header=None, usecols=[0,1])
        # Append the data to the combined DataFrame
        combined_data = pd.concat([combined_data, data], axis=1)

    # Calculate the sum of each column in the combined DataFrame
    total_counts = combined_data.sum(axis=1)

    df = pd.DataFrame({'counts': total_counts,
                                 'percent': [("{a:.2f}%").format(a=total_counts[i]/total_counts['reads']*100) for i in total_counts.index]})
    # Save the total counts to a TSV file
    df.to_csv(output_file, sep='\t', header=None)

def main():
    parser = argparse.ArgumentParser(description="Combine and summarize TSV files.")
    parser.add_argument("file_pattern", help="Pattern to match the TSV files.")
    parser.add_argument("output_file", help="Path to the output TSV file.")
    args = parser.parse_args()

    combine_files(args.file_pattern, args.output_file)

if __name__ == "__main__":
    main()
