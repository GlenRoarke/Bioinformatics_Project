import sys
import os
import re

bamd = sys.argv[1]
output_folder = sys.argv[2]  # New argument for specifying the output folder
bamList = os.listdir(bamd)

output_file = os.path.join(output_folder, "samples.yaml")  # Create the output file path

with open(output_file, "w") as output:  # Open the output file in the specified folder
    output.write("samples:\n")
    for bam in bamList:
        if bam.endswith(".bam"):
            name = re.search(r'(.*).bam', bam).group(1)
            output.write("  " + name + ": " + os.path.join(bamd, bam) + "\n")
