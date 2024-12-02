import argparse
from tp3 import PLAST



#in CLI ,
#enter string of format :
#
#python plast.py -i CGTAGTCGGCTAACCAGCATAACGCTTGTAAACGTAAGAGCCC -db tRNAs.fasta -E 5 -ss 10 -seed '11111111111'


def main():
    # Create an argument parser
    parser = argparse.ArgumentParser(description="PLAST-like tool for sequence searching.")

    # Define arguments
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Input sequence to search.")
    parser.add_argument("-db", "--database", type=str, required=True,
                        help="Database file in FASTA format.")
    parser.add_argument("-E", "--evalue", type=int, default=5,
                        help="E-value threshold for matches (default: 5).")
    parser.add_argument("-ss", "--secondseuil", type=int, default=10,
                        help="Score threshold for high-scoring segments (default: 10).")
    parser.add_argument("-seed", "--seed", type=str, required=True,
                        help="Seed for k-mers (e.g., '11111111111').")

    # Parse the arguments
    args = parser.parse_args()
    #
    print("_____INPUT DATA_____")
    print(f"Input pattern: {args.input}") #
    print(f"Input db file: {args.database}") 
    print(f"Input min E-value: {args.evalue}") #
    print(f"Input ss: {args.secondseuil}")
    print(f"Input seed: {args.seed}") #
    #pattern, seed, match penality, mismatch penality, min e value, ss, db path
    #
    #
    #
    plast = PLAST(args.input, args.seed,5,-4, args.evalue, args.secondseuil, args.database)
    plast.run()


if __name__ == "__main__":
    main()
