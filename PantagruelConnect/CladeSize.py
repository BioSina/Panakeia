"""
@author Sina Beier

sb59@sanger.ac.uk

Copyright is owned by the Sanger Institute (Genome Research Ltd.) and Sina Beier

This is a script for printing tables from cladegraphs.

"""

import argparse
import re


def main():
    """
        Main function
    """
    with open(args.clades) as c:
        with open(args.out, 'w') as o:
            o.write("clade\tcladeSize\tsisterCladeSize\n")
            header = c.readline()
            for line in c:
                line = re.sub("\n", "", line)
                s = re.split("\t", line)
                clade = s[0]
                split2 = re.split(",", s[1])
                split3 = re.split(",", s[2])
                o.write(clade+"\t"+str(len(split2))+"\t"+str(len(split3))+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("clades", help="Clade definition")
    parser.add_argument("out", help="Output table")

    args = parser.parse_args()
    main()