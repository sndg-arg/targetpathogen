
import pickle
import json
import sys
import os.path


def unpickle_and_save_to_json(pickle_file_path):
    with open(pickle_file_path, 'rb') as file:
        unpickled_data = pickle.load(file)
    base_name = os.path.basename(pickle_file_path)  # Get the base name of the file
    name_without_ext = os.path.splitext(base_name)[0]  # Remove the extension
    json_output_path = os.path.join(os.path.dirname(pickle_file_path), f"{name_without_ext}.json")
    with open(json_output_path, 'w') as json_file:
        json.dump(unpickled_data, json_file)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <pickle_file_path>")
    else:
        pickle_file_path = sys.argv[1]
        unpickle_and_save_to_json(pickle_file_path)
