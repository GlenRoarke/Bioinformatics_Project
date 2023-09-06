import argparse
import os
import re

def create_samples_yaml(bam_dir, output_folder):
    bamList = os.listdir(bam_dir)
    output_file = os.path.join(output_folder, "samples.yaml")

    with open(output_file, "w") as output:
        output.write("samples:\n")
        for bam in bamList:
            if bam.endswith(".bam"):
                name = re.search(r'(.*).bam', bam).group(1)
                output.write(f"  {name}: {os.path.join(bam_dir, bam)}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create samples.yaml file for BAM files")
    parser.add_argument("bam_dir", help="Path to the directory containing BAM files")
    parser.add_argument("output_folder", help="Path to the output folder")
    args = parser.parse_args()

    create_samples_yaml(args.bam_dir, args.output_folder)