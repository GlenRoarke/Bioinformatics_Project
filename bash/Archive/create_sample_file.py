#This script will create a "samples.yaml" file to be placed in
#the griffin_GC_and_mappability_correction/config directory
#Usage: python /path/to/recal_bam_files

import sys
import os
import re

bamd = sys.argv[1]
bamList = os.listdir(bamd)

with open("samples.yaml", "w") as output:
	output.write("samples:\n")
	for bam in bamList:
		if bam.endswith(".bam"):
			name = re.search(r'(.*).bam', bam).group(1)
			output.write("  "+name+": "+bamd+"/"+bam+"\n")
