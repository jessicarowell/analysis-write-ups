import os
import sys
import time
from datetime import date
import argparse
import logging as log
from itertools import groupby

def getarg():
    """Set up arguments

    Returns:
        argparse.Namespace: an argparse.Namespace object with all the user-entered arguments
    """
    parser = argparse.ArgumentParser(prog='translate-tr2gen',
                                     usage='%(prog)s [options]',
                                     add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-h', '--help', 
                        action = 'help',
                        help='Show help message and exit\n\n')
    parser.add_argument('-r', '--ref',
                        type=str,
                        required=True,
                        help=('Path to transcript mapping reference file\n'
                        'a tab-delimited file with 4 columns:\n'
                        'transcript name, chromosome name, 0-based start position, CIGAR string\n\n'
                        'e.g., /path/to/mapping_ref.txt\n\n'))
    parser.add_argument('-q', '--query',
                        type=str,
                        required=True,
                        help=('Path to transcript queries\n'
                        'a tab-delimited file with 2 columns: transcript name, 0-based position\n'
                        'e.g., /path/to/queries.txt\n\n'))
    parser.add_argument('-o', '--outputdir',
                        type=str,
                        help=('Output directory\n'
                        'default: current working directory\n\n'),
                        default='./', metavar='')
    parser.add_argument('-f', '--outfile',
                        type=str,
                        help=('Output file name\n'
                        'default: queries_mapped.txt\n\n'),
                        default='queries_mapped.txt', metavar='')
    args = parser.parse_args()
    return args

def setuplog(outdir,logfile):
    """Set up logging

    Args:
        outdir (str): path to output directory 
        logfile (str): path to log file
    """
    with open(logfile, 'w'):
        pass
    handlers = [log.FileHandler(logfile)]
    log.basicConfig(level=log.DEBUG, 
                        format='%(asctime)s - %(message)s',
                        datefmt='%H:%M:%S',
                        handlers=handlers)
    log.info(f'Date: %s', str(date.today()))

def check_path(path):
    """Check if path exists and make sure it ends with "/"

    Args:
        path (str): A directory path

    Returns:
        str: A directory path ending with "/"
    """
    if not os.path.isdir(path):
        try:
            os.mkdir(path)
        except:
            log.error(f'ERROR Directory {path} does not exist and cannot be created. Exiting...')
            exit(1)
    if path[-1] != '/':
        path += '/'
    return path

def checkfile(filename):
    """Check if file path is valid and return it, or exit if not valid.

    Args:
        filename (str): a file path

    Returns:
        str: the file path, or exits if file not found or empty
    """
    if filename:
        if os.path.isfile(filename):
            log.info(f'Found file: {filename}')
            if os.stat(filename).st_size == 0:
                log.error(f'ERROR File: {filename} is empty. Exiting...')
                exit(1)
        else:
            log.error(f'ERROR Cannot locate {filename}. Exiting...')
            exit(1)
    return filename

def read_file(filename):
    """
    Parse a 4-column, tab-delimited file into a dict.
    Store first col as the key, and the last 3 cols as value (as a list).

    Args:
        filename (str): Path of file to parse

    Returns:
        dict: A dictionary where 1st col of the file are keys, last 3 cols are values.
    """
    d = {}
    with open(filename, "r") as fh:
        for line in fh:
            # todo: need some error-checking here - assumes \t and 4 cols
            tr, chrname, rpos, cigar = line.strip().split('\t', 3)
            rpos = int(rpos)
            tr_map_vals = [chrname,rpos,cigar]
            d[tr] = tr_map_vals
    return d


def parse_cigar(cigar, qpos, rpos):
    """
    Given a CIGAR string, a transcript ending coordinate, and a reference starting coordinate, 
    map the transcript to the reference and return the reference sequence ending coordinate.

    Args:
        cigar (str): Cigar string
        qpos (int): Ending coordinate of the query sequence
        rpos (int): Starting coordinate of the reference match

    Returns:
        int: Ending coordinate of the reference match
    """
    q_consume = ("M", "I", "S", "=", "X") # alignment match types that consume query
    ref_consume = ("M", "D", "H", "=", "X") # alignment match types that consume reference
    result = 0 # query position counter (0-based)
    ref_result = rpos  # reference position counter
    cig_iter_grp = groupby(cigar, lambda chr: chr.isdigit()) # Create cigar string groups
    # Loop thru cigar groups and consume query and ref positions
    for _, length_digits in cig_iter_grp:
        length = int(''.join(length_digits)) # alignment match length
        op = next(next(cig_iter_grp)[1]) # alignment match type
        # if we haven't reached the end position, keep processing cigar string
        if result < qpos:
            # we consume the whole match length if it doesn't make query length go past the end position
            if (result + length) <= qpos: 
                # add match length to query/ref seq according to match type
                if op in q_consume:
                    result += length
                if op in ref_consume:
                    ref_result += length
            # if the match length is too long, we consume all/part of it depending on the match type
            else:
                # if match gets consumed by both query and ref, we can only consume up to the end position
                if op in q_consume and op in ref_consume:       
                    leftover = result
                    result += (qpos - leftover)
                    ref_result += (qpos - leftover)
                    return ref_result
                # if match gets consumed only by the ref, we can consume all and continue on
                elif op in ref_consume:
                    ref_result += length               
    return ref_result

if __name__ == '__main__':

    starttime = time.time()
    args = getarg()

    # Define file paths
    # Ensure trailing slash
    outdir = check_path(args.outputdir)

    # Set up log file
    logfile = ''.join([outdir, 'transcript2genomic.log'])
    setuplog(outdir, logfile)
    log.info(f'Command: \n\n\t%s\n', " ".join(sys.argv))
    log.info(f'Arguments: \n\n\t%s\n', str(args.__dict__))

    ref_file = checkfile(args.ref)
    query_file = checkfile(args.query)

    # File 1: Read in, store as a hash table
    try:
        d = read_file(ref_file)
    except:
        log.error(f'ERROR Cannot open and parse {filename}. Please check the format. Exiting...')
        exit(1)

    # Read in file as a list of positions
    try:
        with open(query_file, "r") as fh:
            tr_list = [line.rstrip().split() for line in fh.readlines()]
    except:
        log.error(f'ERROR Cannot parse {query_file}. Please check format. Exiting...')
        exit(1)

    # Convert positions to numeric
    try:
        for i in range(len(tr_list)):
            tr_list[i][1] = int(tr_list[i][1])
    except:
        log.error(f'ERROR Unable to convert position to numeric. Please check column 2 of {query_file}. Exiting...')
        exit(1)

    # Go through list, retrieving genomic coords
    queries = []
    try:
        for tr, tpos in sorted(tr_list, key = lambda i: (i[0], i[1])):
            # find tr in d and return the values (chrname, rpos, cigar)
            chrname, rpos, cigar = d.get(tr)
            # feed to cigar fx the cigar string, tpos, rpos
            mapped_pos = parse_cigar(cigar, tpos, rpos)
            # add (tr, qpos, chr, rpos) to list
            queries.append([tr, tpos, chrname, mapped_pos])
    except:
        # todo: Need some better error-checking here
        log.error(f'ERROR Unable to complete mapping. Exiting...')
        exit(1)

    # Put the queries in same order as original query file
    try:
        keys = {tuple(e): i for i, e in enumerate(tr_list)}
        sorted_queries = sorted(queries, key=lambda x: keys.get(tuple(x[:2])))
    except:
        log.error(f'ERROR Unable to sort list. Exiting...')
        exit(1)

    # Write the list of lists to a file
    with open(outdir + args.outfile, 'w') as fh:
        lines = ['\t'.join(map(str, query)) + '\n' for query in sorted_queries]
        fh.writelines(lines)

    log.info(f'Program finished')
    endtime = time.time()
    log.info(
        f'Execution time (sec): %s', 
        str(round(endtime - starttime))
    )