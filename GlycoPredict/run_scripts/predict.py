from src import glyco_predicter
import time
import argparse

###Run script for using GlycoPredicter in terminal###
def main():
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Run glyco_predicter with arguments.")
    
    # Add arguments as required by glyco_predicter.predict()
    parser.add_argument("--input_path", required=True, help="Path to input file")
    parser.add_argument("--save_name", required=True, help="name to save results under")
    
    # Parse the arguments from the terminal
    args = parser.parse_args()

    start_time = time.time()
    glyco_predicter.predict(args.input_path, args.save_name)
    end_time = time.time()
    print(f"Program finished in: {end_time - start_time} seconds")
