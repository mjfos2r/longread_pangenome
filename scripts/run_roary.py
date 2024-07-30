import os
import glob
import subprocess

# Define the path to the Roary Apptainer image
roary_image = "/home/mf019/longread_pangenome/singularity/images/roary_latest.sif"

# Function to run Roary using Apptainer
def run_roary(input_folder, output_folder):
    gffs = glob.glob(f'{input_folder}/*.gff3')
    command = [
        "apptainer", "exec",
        "--bind", f"{input_folder}:/data/input:r", # bind input to container, give explicit read perms
        "--bind", f"{output_folder}:/data/output:rw", # bind output to container, give explicit rw perms
        roary_image, "roary", # call container image and then run roary inside.
        "-p", "32", "-s", # 32 threads, -s to NOT split paralogs.
        "-e", "--mafft", # prank aln, mafft
        "-f", "/data/output", # Output folder
    ] + [f"/data/input/{os.path.basename(file)}" for file in gffs]
    subprocess.run(command, check=True)

# actually run roary
input_folder = '/home/mf019/longread_pangenome/longread_analysis/paired_assemblies/annotation/longread'
output_folder = '/home/mf019/longread_pangenome/longread_analysis/v5'
# Create the output directory if it doesn't exist
os.makedirs(output_folder, exist_ok=True)
os.chmod(output_folder, 0o777) # ugh change perms so roary can write to output
print(f'{input_folder}\n{output_folder}')
# Run Roary
run_roary(input_folder, output_folder)

print("Roary pangenome analysis completed for all folders.")

#roary -p 32 -s -e --mafft