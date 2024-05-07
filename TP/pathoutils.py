import sys
import pickle
import csv

def main(pickle_file, tsv_file):
    # Load the data from the pickle file
    with open(pickle_file, 'rb') as f:
        data = pickle.load(f)
    
    # Check if the data is a dictionary with lists as values
    if isinstance(data, dict) and all(isinstance(value, list) for value in data.values()):
        # Prepare the header and rows for the TSV file
        header = ['Identifier', 'Compounds']
        
        # Write the data to a TSV file
        with open(tsv_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(header)
            for key, compounds in data.items():
                # Convert the list of compounds to a single string, separated by commas
                compounds_str = ', '.join(compounds)
                writer.writerow([key, compounds_str])
    else:
        print("Error: The data structure is not as expected. Expected a dictionary with lists as values.")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <pickle_file> <tsv_file>")
        sys.exit(1)
    
    pickle_file = sys.argv[1]
    tsv_file = sys.argv[2]
    main(pickle_file, tsv_file)