import sys
import os
import re

def create_sites_yaml(site_dir, output_file):
    fileList = os.listdir(site_dir)

    with open(output_file, "w") as output:
        output.write("site_lists:\n")
        for f in fileList:
            name = re.search(r'(.*).txt', f).group(1)
            path = os.path.join(site_dir, f)
            output.write(f"  {name}: {path}\n")

if __name__ == "__main__":
    site_dir = sys.argv[1]
    output_file = sys.argv[2]

    create_sites_yaml(site_dir, output_file)
