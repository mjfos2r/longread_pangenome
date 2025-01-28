import os
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

def get_input_files(input_dir):
    # Get a sorted list of files in the input directory
    return sorted([os.path.join(input_dir, f) for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))])

def set_output(file1, file2, output_dir):
    # Create an output path based on the two input files
    output_filename = f"{os.path.basename(file1)}_vs_{os.path.basename(file2)}.tsv"
    return os.path.join(output_dir, output_filename)

def get_alignment_command(file1, file2, output_path, prog):
    # Construct the alignment command based on the chosen program
    if prog == "mauve":
        return get_mauve_command(file1, file2, output_path)
    elif prog == "mummer":
        return get_mummer_command(file1, file2, output_path)

def run_alignment(command, prog):
    # Execute the alignment command
    if prog == "mauve":
        run_mauve(command)
    elif prog == "mummer":
        run_pgv_mummer(command)
    return f"Finished alignment: {command}"

def get_mauve_command(file1, file2, output_path):
    # Placeholder for generating a Mauve command
    return ["mauveAligner", "--output", output_path, "--ref", file1, file2]

def run_mauve(command):
    # Placeholder for running a Mauve command
    os.system(' '.join(command))

def get_mummer_command(file1, file2, output_path):
    # Placeholder for generating a Mummer command
    return ["pgv-", "--prefix", output_path, file1, file2]

def run_pgv_mummer(command):
    # Placeholder for running a Mummer command
    os.system(' '.join(command))

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="A script to run alignments on a directory of assemblies in an all-vs-all fashion.")
    # Add the arguments
    parser.add_argument("prog", type=str, help='type of alignment to run [Mauve or Mummer]')
    parser.add_argument('input_dir', type=str, help='The directory of inputs')
    parser.add_argument('output_dir', type=str, help='The directory for outputs')
    parser.add_argument("--max_workers", type=int, default=None, help="Maximum number of workers to run in parallel")
    # Parse the arguments
    args = parser.parse_args()
    
    # Print the arguments
    print(f"alignment type {args.prog}")
    print(f"Input dir {args.input_dir}")
    print(f"Output dir: {args.output_dir}")

    input_files = get_input_files(args.input_dir)

    with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        futures = []
        for i, file1 in enumerate(input_files):
            for j, file2 in enumerate(input_files):
                if i < j:  # Ensure each pair is only processed once and exclude self-comparisons
                    output_path = set_output(file1, file2, args.output_dir)
                    command = get_alignment_command(file1, file2, output_path, args.prog)
                    futures.append(executor.submit(run_alignment, command, args.prog))
        
        with tqdm(total=len(futures)) as pbar:
            for future in as_completed(futures):
                try:
                    result = future.result()
                    tqdm.write(result)  # Use tqdm.write to print stdout while keeping the progress bar
                except Exception as e:
                    tqdm.write(f"Error: {e}")
                pbar.update(1)

if __name__ == "__main__":
    main()