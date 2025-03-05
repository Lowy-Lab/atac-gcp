import argparse
import sys

def main():
    if sys.argv[1] == "human":
        chroms={"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"}
    elif sys.argv[1] == "mouse":
        chroms = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY"}
    else:
        print("Currently only mouse and human are supported")
        sys.exit(1)
    chrom_list = []
    for n in sys.stdin.read().split(" "):
        if n not in chroms:
            chrom_list.append(n)
    print(" ".join(chrom_list))
     
     


if __name__ == '__main__':
  main()