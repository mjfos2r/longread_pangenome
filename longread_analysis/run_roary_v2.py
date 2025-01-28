import os
import glob
import subprocess

# Get CWD
cwd = os.getcwd()

# Define the path to the Roary Apptainer image
roary_image = "/home/mf019/longread_pangenome/singularity/images/roary_latest.sif"

# Function to run Roary using Apptainer
def run_roary(input_folder, output_folder):
    gffs = glob.glob(f'{input_folder}/*.gff3')
    if not gffs:
        raise FileNotFoundError("No GFF3 files found in the input folder.")
    
    print("GFF3 files to be processed:", gffs)
    
    # Check permissions inside the container
    check_command = [
        "apptainer", "exec",
        "--bind", f"{input_folder}:/data/input:rw",
        "--bind", f"{output_folder}:/data/output:rw",
        roary_image, "sh", "-c",
        "ls -ld /data/input /data/output && touch /data/output/testfile"
    ]
    
    try:
        subprocess.run(check_command, check=True)
        print("Directory permissions and writability check passed inside the container.")
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        print("Command output:", e.output)
        raise
    
    # Roary command
    command = [
        "apptainer", "exec",
        #"--bind", f"{input_folder}:/data/input:rw", # bind input to container, give explicit read perms
        #"--bind", f"{output_folder}:/data/output:rw", # bind output to container, give explicit rw perms
        roary_image, "roary", # call container image and then run roary inside.
        "-p", "32", "-s", # 32 threads, -s to NOT split paralogs.
        "-e", "--mafft", # prank aln, mafft
        "-f", output_folder, #"/data/output", # Output folder
        ] + [f"{input_folder}/{os.path.basename(file)}" for file in gffs]
    
    print("Running command:", " ".join(command))
    
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        print("Command output:", e.output)
        raise

# Actually run Roary
input_folder = '/home/mf019/longread_pangenome/longread_analysis/paired_assemblies/annotation/longread'
output_folder = '/home/mf019/longread_pangenome/longread_analysis/v5'

# Create the output directory if it doesn't exist
os.makedirs(output_folder, exist_ok=True)
os.chmod(output_folder, 0o777) # Explicitly change perms so roary can write to output

# Check permissions of input and output directories
input_folder_perms = os.stat(input_folder).st_mode
output_folder_perms = os.stat(output_folder).st_mode

print(f'Input folder: {input_folder} (permissions: {oct(input_folder_perms)})')
print(f'Output folder: {output_folder} (permissions: {oct(output_folder_perms)})')

# Run Roary
run_roary(input_folder, output_folder)

print("Roary pangenome analysis completed for all folders.")

