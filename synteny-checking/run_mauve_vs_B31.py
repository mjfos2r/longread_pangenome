import os
import glob
import argparse
import subprocess

def get_inputs(in_dir):
    genomes = glob.glob(f'{in_dir}/*.fasta')
    return genomes

def set_output(in_file, out_dir):
    isolate_name = os.path.basename(in_file).split('.')[0]
    os.makedirs(f'{out_dir}/{isolate_name}', exist_ok=True)
    output_file = f'{out_dir}/{isolate_name}/{isolate_name}_vs_B31'
    return output_file

def get_command(in_file, out_file, ref_genome_path):
    # Set up the command to run one genome vs the B31 reference.
    progressive_mauve_path = '/Applications/Mauve.app/Contents/MacOS/progressiveMauve'
    command = [f'{progressive_mauve_path}', f'--output={out_file}.xmfa',
                 f'--backbone-output={out_file}.backbone',
                 f'--output-guide-tree={out_file}.gtree',
                 f'{ref_genome_path}', f'{in_file}']
    return command

def run_mauve(command):
    # Run the command
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print("Mauve alignment completed successfully.")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print("An error occurred while running Mauve.")
        print(e.stderr)

def parse_results(results):
    # Parse out the useful results?
    # visualize?
    pass

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="A script to run mauve on a directory of assemblies against the B31 reference genome.")
    # Add the arguments
    parser.add_argument('input_dir', type=str, help='The directory of inputs')
    parser.add_argument('output_dir', type=str, help='The directory for outputs')
    parser.add_argument("reference", type=str, help="path to reference genome for alignment")
    #parser.add_argument('cores_for_mauve', type=int, help='how many cores to run mauve with')
    # Parse the arguments
    args = parser.parse_args()

    # Print the arguments (for demonstration purposes)
    print(f"Input dir {args.input_dir}")
    print(f"Output dir: {args.output_dir}")
    print(f"reference: {args.reference}")
    #print(f"Cores for Mauve: {args.cores_for_mauve}")

    genomes = get_inputs(args.input_dir)
    for genome in genomes:
        out_file = set_output(genome, args.output_dir)
        command = get_command(genome, out_file, args.reference)
        run_mauve(command)
        print(f"Mauve against B31 for {genome} finished!")

if __name__ == "__main__":
    main()
