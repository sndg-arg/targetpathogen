import argparse
import pandas as pd
import os

def main(csv_file):
    # Read the CSV file using pandas
    data = pd.read_csv(csv_file)
    data["locus_tag"] = [str(x).split(" ")[-1] for x in data["Gene Names"]]
    data = data[["Entry", "locus_tag"]]
    data = data.dropna(subset=["Entry"])
    data = data.dropna(subset=["locus_tag"]) 
    # Do further processing or analysis with the data
    print("Data loaded successfully:\n", data)
    # Save the data as csv
    directory, filename = os.path.split(csv_file)
    output_file = os.path.join(directory, "uniprot_locustag.csv")
    data.to_csv(output_file, index=False)

if __name__ == "__main__":
    # Creating ArgumentParser object
    parser = argparse.ArgumentParser(description="Load a CSV file using pandas")
    
    # Adding a positional argument for the CSV file
    parser.add_argument("csv_file", help="Path to the CSV file")

    # Parsing arguments
    args = parser.parse_args()

    # Call the main function with the provided CSV file path
    main(args.csv_file)