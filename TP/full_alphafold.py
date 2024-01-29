import argparse
import csv
import subprocess

def run_alphafold(uniprot_entry, locus_tag):
    command = [
        "python", "-m", "TP.alphafold",
        "-pr", "opt/p2rank/distro/prank",
        "-o", "data/025/NC_002516.2/NC_002516.2/alphafold",
        "-T", "10",
        "-ltag", locus_tag,
        "-nc",
        "-c", uniprot_entry
    ]
    subprocess.run(command)

def main(csv_file):
    with open(csv_file, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header
        for row in reader:
            uniprot_entry, locus_tag = row
            run_alphafold(uniprot_entry, locus_tag)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run AlphaFold with UniProt entry codes and locus tags from a CSV file")
    parser.add_argument("csv_file", help="Path to the CSV file containing UniProt entry codes and locus tags")
    args = parser.parse_args()
    main(args.csv_file)