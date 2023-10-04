# import numpy library for flexible arrays
import numpy as np

#nested loops are bad. Use methods to break out the function.
# we use two here - the first pulls out the coordinates for a fragment
# the second increments the current array

def increment_alignment(start, end, current_array):
    #In python, better to ask for forgiveness than permission, so encode this in a
    #try/except block to handle out of range values
    try:
        for i in range(start, end+1):
            current_array[i]+=1
    except IndexError as e:
        print(e)

def process_line(line, current_LRR, current_array):
    alignment_attributes = line.split()
    LRR_name = alignment_attributes[0]

    #set up a new array and a new name
    if current_LRR != LRR_name:
        LRR_length = int(alignment_attributes[1])  # function int converts strings into integers
        current_array = np.zeros(LRR_length)  # set up the array of zeros to match LLR contig length
        current_LRR = LRR_name

    # don't need the 'else' because you want to do the same thing in both cases, so we can move that outside the if statement
    alignment_start = int(alignment_attributes[2])
    alignment_end = int(alignment_attributes[3])
    increment_alignment(alignment_start, alignment_end, current_array)
    return current_LRR, current_array


def get_regions_of_zero(arr):
    """
    Identifies regions in an array where the value is 0.
    :param arr: A numpy array of numbers
    :return: Returns a list of sets of two indices. The first marks where a region of zeros starts, the second marks where the region of zeros stops.
    """

    #Finds ALL indices where the value is 0 (including adjacent indices).
    indices = np.where(arr == 0)[0]

    #if there aren't any indices with a value of 0 in the array, return an empty list
    if len(indices) == 0:
        return []

    #process the list of indices to find the first and the last index for a region of adjacent 0s.
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

all_LRR_breakages = []

def main():
    # set current LRR to empty
    current_LRR = ''
    current_array = None
    LRR_breakages = ''
    all_LRR_breakages.append(LRR_breakages)
    
    # open alignment file for single long-read cluster rep; read each line and assign variables to each part
    with open("LRR_1.txt", "r") as handle:
        for line in handle.readlines():
            line = line.strip()
            #ignore blank lines
            if line:
                current_LRR, current_array = process_line(line, current_LRR, current_array)
                regions_of_zero = get_regions_of_zero(current_array)
                #print(current_array)
                #print(regions_of_zero)
                LRR_breakages = current_LRR, regions_of_zero
                
with open('output.csv', 'w') as handle:
    for LRR_breakages in all_LRR_breakages:
        handle.write(f'{LRR_breakages}\n')
                #all_LRR_breakages.append(LRR_breakages)
                #print(LRR_breakages)

            #write the output to a csv file:
        
    #with open('output.csv', 'w') as handle:
        #handle.write('LRR_name, breakages\n'):
    #    handle.write(f'{all_LRR_breakages}\n')

if __name__ == "__main__":
    main()