import argparse
import numpy as np
import matplotlib.pyplot as plt



def main(args):
    matrix = np.load(args.npy)
    plt.imshow(matrix, cmap='viridis', aspect='auto')
    plt.colorbar(label='Values')
    plt.savefig(args.output, dpi=300, bbox_inches='tight', format='pdf')




if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-i', '--npy', type=str,
                        help="File containing pairwise distance matrix",
                        required=True)

    parser.add_argument('-o', '--output', type=str,
                        help='File containing output result',
                        required=True)
    main(parser.parse_args())