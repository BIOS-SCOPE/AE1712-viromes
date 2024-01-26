# import numpy library for flexible arrays
import numpy as np
import pandas as pd
import argparse
import logging
import pathlib


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--alignment_file", type=str, required=True)
    parser.add_argument("--out_file", type=str, required=True)
    parser.add_argument("--coverage_file", type=str, default='coverage.txt')
    parser.add_argument("--min_pct_align", type=float, default=70)
    # Parse arguments
    return parser.parse_args()


# nested loops are bad. Use methods to break out the function.
# we use two here - the first pulls out the coordinates for a fragment
# the second increments the current array

def increment_alignment(start, end, current_array):
    # In python, better to ask for forgiveness than permission, so encode this in a
    # try/except block to handle out of range values
    try:
        for i in range(start -1, end):
            current_array[i] += 1
    except IndexError as e:
        logging.error(e)


def parse_regions_of_zero_into_rows(regions_of_zero, current_LRR, results):
    for region in regions_of_zero:
        row_list = [current_LRR, region[0], region[1], region[1] - region[0]+1]
        results.append(row_list)


def process_group(name, group, results,coverage_file):
    target_id, target_length = name
    current_array = np.zeros(target_length)  # set up the array of zeros to match LLR contig length

    for i, row in group.iterrows():
        alignment_length = row['target_end'] - row['target_start']
        pct_alignment = alignment_length / row['query_length'] * 100
        if pct_alignment >= args.min_pct_align:
            logging.info(
                f"Found a short read contig that aligns to {target_id} in the region {row['target_start']}-{row['target_end'] }, so incrementing the coverage in this region")
            increment_alignment(row['target_start'], row['target_end'], current_array)
        else:
            logging.debug(f"Short read contig {row['query_id']} aligned to {target_id} in the region {row['target_start']}-{row['target_end']}, but did not meet minimum requirements ({pct_alignment:.1f}%")

    with open(coverage_file, 'a') as handle:
        handle.write(f'{target_id}: {"".join([str(int(x)) for x in current_array.tolist()])}\n')
    regions_of_zero = get_regions_of_zero(current_array)
    parse_regions_of_zero_into_rows(regions_of_zero, target_id, results)


def get_regions_of_zero(arr):
    """
    Identifies regions in an array where the value is 0.
    :param arr: A numpy array of numbers
    :return: Returns a list of sets of two indices. The first marks where a region of zeros starts, the second marks where the region of zeros stops.
    """

    # Finds ALL indices where the value is 0 (including adjacent indices).
    indices = np.where(arr == 0)[0]

    # if there aren't any indices with a value of 0 in the array, return an empty list
    if len(indices) == 0:
        return []

    # process the list of indices to find the first and the last index for a region of adjacent 0s.
    ranges = []
    start = end = indices[0]
    for i in indices[1:]:
        if i - end == 1:
            end = i
        else:
            ranges.append((start, end))
            start = end = i
    ranges.append((start, end))
    return ranges


def main(args):
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

    ## delete the existing coverage file if it exists
    pathlib.Path(args.coverage_file).unlink(missing_ok=True)


    results = []

    # parse file
    minimap2_out = pd.read_csv(args.alignment_file, sep='\t',
                               header=None, usecols=[0, 1, 2, 3, 5, 6, 7, 8, 9, 10],
                               names=['query_id', 'query_length', 'query_start', 'query_end', 'target_id', 'target_length','target_start', 'target_end', 'num_residue_match', 'alignment_block_length'])

    contig_groups = minimap2_out.groupby(['target_id', 'target_length'])
    for name, group in contig_groups:
        process_group(name, group, results, args.coverage_file)

    # now we use pandas magic
    df = pd.DataFrame(results, columns=['target_id', 'start', 'end', 'breakage_length']).merge(minimap2_out[['target_id', 'target_length']].drop_duplicates())
    df.to_csv(args.out_file, sep='\t', index=False)



if __name__ == "__main__":
    args = parse_args()
    main(args)
