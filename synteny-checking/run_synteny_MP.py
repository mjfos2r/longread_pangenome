import os
import glob
import argparse
import subprocess
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

def get_fasta_inputs(in_dir):
    genomes = glob.glob(f'{in_dir}/*.fasta')
    return genomes

def get_gb_inputs(in_dir):
    genomes = glob.glob(f'{in_dir}/*.gbff')
    return genomes

def set_output(in_file, out_dir):
    isolate_name = in_file.split('/')[-1].replace('.gbff','')
    os.makedirs(f'{out_dir}/{isolate_name}', exist_ok=True)
    output_path = f'{out_dir}/{isolate_name}'
    return output_path

def set_multi_output(file1, file2, out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    plasmid1_id = os.path.basename(file1).replace('.gbff','')
    plasmid2_id = os.path.basename(file2).replace('.gbff','')
    filename = f'{plasmid1_id}_vs_{plasmid2_id}'
    #print(output_path)
    return os.path.join(out_dir, filename)

def get_alignment_command(file1, file2, output_path, prog):
    # Construct the alignment command based on the chosen program
    if prog == "mauve":
        return get_mauve_command(file1, file2, output_path)
    elif prog == "mummer":
        return get_mummer_command(file1, file2, output_path)

def get_mauve_command(file1, file2, output_path):
    # Set up the command to run one genome vs the B31 reference.
    isolate_name = file1.split('/')[-1].replace('.gbff','')
    out_file = f'{output_path}/{isolate_name}'
    progressive_mauve_path = "progressiveMauve"
    command = [f'{progressive_mauve_path}', f'--output={out_file}/{isolate_name}.xmfa',
                 f'--backbone-output={out_file}.backbone',
                 f'--output-guide-tree={out_file}.gtree',
                 f'{file2}', f'{file1}']
    return command

def get_mummer_command(file1, file2, output_path):
    # Set up the command to run one genome vs the B31 reference.
    command = ['pgv-mummer', '-q', '--threads=4', f'{file1}', f'{file2}', f'-o {output_path}',
               '--formats=html', '--length_thr=250', '--identity_thr=75', '--show_scale_bar', '--curve',]
    return command

def run_mauve(command):
    # Run the command
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    return result.stdout if result.returncode == 0 else result.stderr

def run_pgv_mummer(command):
    # Run the command
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    #print("pgv-mummer alignment completed successfully.")
    return result.stdout if result.returncode == 0 else result.stderr

def parse_results(results):
    # Parse out the useful results?
    # visualize? ## ITS LITERALLY ALREADY DOING THIS :)
    pass

#def run_alignment(prog, genome, output_path, reference):
#    if prog == "mauve":
#        command = get_mauve_command(genome, output_path, reference)
#        run_mauve(command)
#        return f"progMauve against {os.path.basename(reference)} for {genome} finished!"
#    elif prog == "mummer":
#        command = get_mummer_command(genome, output_path, reference)
#        run_pgv_mummer(command)os.path.basename(
#   )gv-mummer command for {genome} against {os.path.basename(reference)} finished!"
def run_alignment(command, prog):
    # Execute the alignment command
    try:
        if prog == "mauve":
            run_mauve(command)
        elif prog == "mummer":
            run_pgv_mummer(command)
        return f"Finished alignment: {' '.join(command)}"
    except Exception as e:
        return f"Error running {' '.join(command)}\nException:{e}"

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="A script to run mauve on a directory of assemblies against the B31 reference genome.")
    # Add the arguments
    parser.add_argument("prog", type=str, help='type of alignment to run [Mauve or Mummer]')
    parser.add_argument('input_dir', type=str, help='The directory of inputs')
    parser.add_argument('output_dir', type=str, help='The directory for outputs')
    #parser.add_argument("reference", type=str, help="path to reference genome for alignment")
    #parser.add_argument('cores_for_mauve', type=int, help='how many cores to run mauve with')
    # Parse the arguments
    args = parser.parse_args()

    # Print the arguments
    print(f"alignment type {args.prog}")
    print(f"Input dir {args.input_dir}")
    print(f"Output dir: {args.output_dir}")
    #print(f"reference: {args.reference}")
    #ref_name = args.reference.split('/')[-1].replace('.gbff','')

    if args.prog == "mauve":
        genomes = get_fasta_inputs(args.input_dir)
    elif args.prog == "mummer":
        if args.input_dir.endswith('.gbff'):
            genomes = [args.input_dir]
        else:
            genomes = get_gb_inputs(args.input_dir)
    else:
        raise ValueError("Unsupported alignment type. Use 'Mauve' or 'Mummer'.")

    with ProcessPoolExecutor(max_workers=90) as executor: # max of 11 workers at 8 cores a pop :)
        futures = []
        for i, genome1 in enumerate(genomes):
            for j, genome2 in enumerate(genomes):
                if i < j: # process each pair only once, do not compare against self.
                    output_path = set_multi_output(genome1, genome2, args.output_dir)
                    command = get_alignment_command(genome1, genome2, output_path, args.prog)
                    futures.append(executor.submit(run_alignment, command, args.prog))
            # Take the list of futures, add an alignment to it, give it 1. the function to call, and then the args for that function.
        with tqdm(total=len(futures)) as pbar:
            for future in as_completed(futures):
                try:
                    result = future.result()
                    tqdm.write(result)
                except Exception as e:
                    tqdm.write(f"Error: {e}")
                pbar.update(1)
    ##print(f"Cores for Mauve: {args.cores_for_mauve}")
    #match args.prog:
    #    case "mauve":
    #        genomes = get_fasta_inputs(args.input_dir)
    #        for genome in genomes:
    #            output_path = set_output(genome, args.output_dir)
    #            command = get_mauve_command(genome, output_path, args.reference)
    #            run_mauve(command)
    #            print(f"progMauve against {ref_name} for {genome} finished!")
    #    case "mummer":
    #        if args.input_dir.endswith('.gbff'):
    #            genomes = args.input_dir
    #        else:
    #            genomes = get_gb_inputs(args.input_dir)
    #        for genome in genomes:
    #            output_path = set_output(genome, args.output_dir)
    #            command = get_mummer_command(genome, output_path, args.reference)
    #            print(' '.join(command))
    #            print(f"running pgv-mummer for {genome} against {ref_name}!")
    #            run_pgv_mummer(command)
    #            print(f"pgv-mummer command finished!\n")
#
if __name__ == "__main__":
    main()
