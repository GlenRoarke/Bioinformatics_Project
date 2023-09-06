#This script generates a sites.yaml file from a given directory path
#Example usage: python create_sites_file.py /user/work/ec20449/Project24/Griffin/Ref/30000_unfiltered_sites_CIS_BP_v2

import sys
import os
import re

siteDir = sys.argv[1]
fileList = os.listdir(siteDir)

output = open("sites.yaml", "w")
output.write("site_lists:\n")
for f in fileList:
    name = re.search(r'(.*).txt', f).group(1)
    path = siteDir + "/" + f
    output.write("  "+name+": "+path+"\n")
output.close() 



