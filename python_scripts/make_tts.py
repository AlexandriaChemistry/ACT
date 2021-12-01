#!/usr/bin/env python3

"""Make train and tests sets from the Alexandria Library.

With 2 numbers being given by the user, this script generates a .dat file
with train and test compounds.
By default, it will place the sets in AlexandriaLib/traintestsets. The
directory will be created if it does not exist.
"""

import os
import sys
import argparse as ap
from random import shuffle


def parse_arguments() -> ap.Namespace:
    """Parses command-line arguments."""
    parser: ap.ArgumentParser = ap.ArgumentParser()

    parser.add_argument('-ntr', '--ntrain', type=int, default='',
                        help='Number of compounds in training set',
                        required=True)
    parser.add_argument('-nte', '--ntest', type=int, default='',
                        help='Number of compounds in test set',
                        required=True)
    parser.add_argument('-alp', '--alib_path', type=str, default='',
                        help='Directory where the Alexandria Library is located. Default: $AlexandriaLib')
    parser.add_argument('-dest', '--destination', type=str, default='',
                        help='Directory to store the set. Default: $AlexandriaLib/traintestsets/')
    parser.add_argument('-alpha', '--alphabetical', action='store_true',
                        help='Sort compounds alphabetically in the output file')              
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print information in console')

    args: ap.Namespace = parser.parse_args()
    return args


def main() -> None:
    """Execute main flow!"""
    args: ap.Namespace = parse_arguments()
    
    # Check validity of the arguments
    if args.ntrain <= 0:
        msg = '\nThe number of training samples must be positive!'
        sys.exit(msg)
    if args.ntest < 0:
        msg = '\nThe number of test samples must be nonnegative!'
        sys.exit(msg)

    # If verbose, print provided arguments
    if args.verbose:
        print(f"\n{'Number of TRAIN compounds:' : <27} {args.ntrain : <35}")
        print(f"{'Number of TEST compounds:' : <27} {args.ntest : <35}")
        print(f"{'Provided AlexandriaLib:': <27} {args.alib_path : <35}")
        print(f"{'Provided destination:': <27} {args.destination : <35}")

    # Check if environmental variable AlexandriaLib or the provided Library path exist
    env_var = 'AlexandriaLib'
    if args.alib_path:
        if not os.path.exists(args.alib_path):
            msg: str = f'\nThe path {args.alib_path} does not exist!'
            sys.exit(msg)
        if not os.path.isdir(args.alib_path):
            msg: str = f'\nThe path {args.alib_path} is not a directory!'
            sys.exit(msg)
    else:
        if env_var not in os.environ:
            msg: str = f'\n<{env_var}> environmental variable is not defined!'
            msg += '\nPlease define it and run this script again...'
            sys.exit(msg)

    # Get path for the environmental variable
    alib_path = os.getenv(env_var) if not args.alib_path else args.alib_path
    if args.verbose:
        msg: str = f'\nAlexandria Library should be located at'
        msg += f' {alib_path}'
        print(msg)

    # Check if path to the Library exists (will exist if provided by the user
    # and made it this far).
    # If it exists, move to it. Otherwise, halt!
    if not os.path.exists(alib_path):
        msg: str = f'\nThe path {alib_path} does not exist!'
        sys.exit(msg)
    if not os.path.isdir(alib_path):
        msg: str = f'\nThe path {alib_path} is not a directory!'
        sys.exit(msg)
    if args.verbose:
        msg: str = f'\nChanging working directory to {alib_path}...'
        print(msg)
    os.chdir(alib_path)

    # Check that there exists a <compounds> directory in the lib
    compounds_path = 'compounds'
    if not os.path.exists(compounds_path):
        msg: str = f'\nThe path {alib_path}/{compounds_path} does not exist!'
        sys.exit(msg)
    if not os.path.isdir(compounds_path):
        msg: str = f'\nThe path {alib_path}/{compounds_path} is not a directory!'
        sys.exit(msg)

    # Check if there are any compounds in the <compounds directory>
    compounds: list = os.listdir(compounds_path)
    n_compounds: int = len(compounds)
    if not compounds:
        msg: str = f'\nThere are no compounds in {alib_path}/{compounds_path}!'
        sys.exit(msg)
    if args.verbose:
        msg: str = f'\nThere are {n_compounds} compounds in the library.'
        print(msg)

    # Check if there are enough compounds to satisfy the needs of the user
    n_train_and_test: int = args.ntrain + args.ntest
    if n_train_and_test > n_compounds:
        msg: str = f'\nYou requested ntrain {args.ntrain} and ntest {args.ntest}'
        msg += f' (total {n_train_and_test}), but'
        msg += f' the library only has {n_compounds} compounds!'
        sys.exit(msg)

    # Shuffle the list of compounds
    if args.verbose:
        msg: str = '\nShuffling the compounds...'
        print(msg)
    shuffle(compounds)

    # Assign training and testing labels
    ttset: list = [f'{compounds[i]}|Train\n' for i in range(args.ntrain)]
    ttset += [f'{compounds[i]}|Test\n' for i in range(args.ntrain, n_train_and_test)]

    # Sort alphabetically if needed
    if args.alphabetical:
        if args.verbose:
            msg: str = '\nSorting compounds alphabetically...'
            print(msg)
        ttset.sort(key=lambda name: name.split('|')[0])
        

    # Get desination path to write and move into it
    ttsets_def_path: str = 'traintestsets'
    if args.destination:
        if not os.path.exists(args.destination):
            msg: str = f'\nThe destination path {args.destination} does not exist!'
            sys.exit(msg)
        if not os.path.isdir(args.destination):
            msg: str = f'\nThe destination path {args.destination} is not a directory!'
            sys.exit(msg)
    else:
        # Create the <traintestsets> directory if it does not exist
        if not os.path.exists(ttsets_def_path):
            if args.verbose:
                msg: str = '\nCreating <traintestsets> directory...'
                print(msg)
            os.mkdir(ttsets_def_path)
    
    ttsets_path: str = f'{alib_path}/{ttsets_def_path}' if not args.destination else args.destination
    if args.verbose:
        msg: str = f'\nChanging working diretory to {alib_path}/{ttsets_def_path}...'
        print(msg)
    os.chdir(ttsets_path)

    # Get existing sets
    existing_sets: list = os.listdir('.')
    # How many datasets with the amount of ntrain and ntest already exist?
    existing_versions: list = [int(eset.split('_')[-1].split('.')[0][1:])
                               for eset in existing_sets
                               if (int(eset.split('_')[1]) == args.ntrain
                                   and int(eset.split('_')[3]) == args.ntest)]
    last_version: int = 0 if not existing_versions else max(existing_versions)
    
    # Create file name and write!
    fname: str = f'ntrain_{args.ntrain}_ntest_{args.ntest}_v{last_version+1}.dat'
    if args.verbose:
        msg: str = f'\nWriting {fname}...'
        print(msg)
    with open(fname, 'w') as f:
        f.writelines(ttset)
    
    if args.verbose:
        msg: str = '\nDONE!\n'
        print(msg)


if __name__ == '__main__':
    main()