"""
Author: Matthew Harper
File: ExpanGe.py

The main python script for the ExpanGe program.
ExpanGe is a bioinformatics tool that is used to analyze a comparison of a series of mums
and identify whether there is an expansion or contraction of the distances between the mums
and whether an inversion occurs in the sequence.
"""

import sys
import getopt
import pandas as pd
from util import Gene, GeneMap

"""
TO DO:
- Implement way to detect head and tail of inversions for proper calculations (DONE)
- Add way to splice away the header and add on new columns (DONE)
- Clean up code (add comments, remove unnecessary code, remove old to do statements and comments)
- Test program with larger file to determine runtime
- Write display_help information
- Write the Readme for the project
- See if way to clean up the look of the output file
- Find way to deal with outlier of data
    - Diagonal chaining
    - Average distance
- Fix inverted distance calculations (DONE)
- Gene Map implementation (EMAIL)
- Make universal way to detect transpositions (DONE)
"""


def identify_inversions(gene_sequence):
    """
    identify_inversions is a function designed to read through the lines in the .coord file and process all the data for
    later use. It identifies where the reversals occur and will set the values for any values that are to be ignored.

    :param gene_sequence: The list holding all the values from the .coords file
    :return: None
    """
    inversion_count = 0
    inv_flag = False
    for x in range(len(gene_sequence)):
        if x != 0:
            current = gene_sequence[x]
            if current.ignore is False:

                ignore_flag = False
                end_flag = False

                if x == (len(gene_sequence) - 1):
                    return None

                iter = 1
                prev = gene_sequence[x-iter]

                while ignore_flag is False and end_flag is False:
                    if prev.ignore is False:
                        ignore_flag = True
                    else:
                        iter = iter + 1
                        if (x-iter) >= 0:
                            end_flag = True
                        else:
                            prev = gene_sequence[x-iter]

                if end_flag is False:
                    start2 = prev.start2
                    start1 = current.start2

                    if start2 > start1:
                        if inv_flag is False:
                            inversion_count = inversion_count + 1
                            inv_flag = True
                            head_index = find_prev_valid(gene_sequence, x)
                            gene_sequence[head_index].inv_head = True
                            gene_sequence[head_index].reversed = True
                        gene_sequence[x].reversed = True
                        gene_sequence[x].inv_count = inversion_count
                    else:
                        if inv_flag is True:
                            gene_sequence[x].inv_tail = True
                            last_index = find_prev_valid(gene_sequence, x)
                            gene_sequence[last_index].inv_last = True
                        inv_flag = False


def find_prev_valid(genes, index):
    """
    Search function that will take in a node and find the closest, previous mum that is valid.
    :param genes: The gene sequence
    :param index: Index of current mum
    :return: The previous, valid mum
    """
    new_index = index - 1

    if new_index < 0:
        return None

    while new_index >= 0:
        if genes[new_index].ignore is False:
            return new_index
        new_index = new_index - 1

    return None


def find_next_not_inverted(genes, index):
    """
    Search function that will find the next mum in the data that is not inverted
    :param genes: The gene sequence
    :param index: The index of the current mum
    :return: the index of the next not inverted mum
    """
    new_index = index + 1

    if new_index == len(genes):
        return None

    while new_index < len(genes):
        if genes[new_index].ignore is False and genes[new_index].reversed is False:
            return new_index
        new_index = new_index + 1

    return None


def find_prev_not_inverted(genes, index):
    """
    Search function that will find the previous mum that is not inverted.
    :param genes: The gene sequence
    :param index: The index of the current mum
    :return: The index of the prev not inverted mum
    """
    new_index = index - 1

    if new_index < 0:
        return None

    while new_index >= 0:
        if genes[new_index].ignore is False and genes[new_index].reversed is False:
            return new_index
        new_index = new_index - 1

    return None


def find_next_valid(genes, index):
    """
    Search function that finds the next valid node in the linked list
    :param genes: Gene sequence
    :param index: Current index
    :return: Index of next valid mum
    """
    new_index = index + 1

    if new_index == len(genes):
        return None

    while new_index < len(genes):
        if genes[new_index].ignore is False:
            return new_index
        new_index = new_index + 1

    return None


def cutoff_check(val1, val2, cutoff):
    """
    Cutoff check function will check if the absolute value of delta_r or delta_q is
    greater than the inputted cutoff value. This is to limit the number of same chromosome
    tranpositions that are taken into account for the calculations
    """
    if abs(val1) >= cutoff or abs(val2) >= cutoff:
        return True

    return False


def calculate_distances(gene_sequence, cutoff=1e6):
    """
    Function to calculate the distances between each gene in the coords file.
    If a gene's previous gene is invalid, it will find the closest valid gene to calculate the
    distance from.
    Utilizes a recursive helper function.
    :param gene_sequence: The data from the coords file
    :param cutoff: Cutoff value
    :return: None
    """

    for x in range(len(gene_sequence)):
        if gene_sequence[x].ignore is False:
            # Distance calculation for non inverted mum
            if gene_sequence[x].reversed is False:
                # Distance calculation for tail of the inversion
                if gene_sequence[x].inv_tail is True:
                    non_inv = find_prev_not_inverted(gene_sequence, x)
                    prev_index = find_prev_valid(gene_sequence, x)
                    if prev_index is not None and non_inv is not None:
                        gene_sequence[x].delta_r = gene_sequence[x].start1 - gene_sequence[prev_index].end1
                        gene_sequence[x].delta_q = gene_sequence[x].start2 - gene_sequence[non_inv].end2
                        gene_sequence[x].delta_x = gene_sequence[x].delta_q - gene_sequence[x].delta_r

                        flag = cutoff_check(gene_sequence[x].delta_r, gene_sequence[x].delta_q, cutoff)
                        if flag:
                            gene_sequence[x].delta_r = "--"
                            gene_sequence[x].delta_q = "--"
                            gene_sequence[x].delta_x = "--"

                else:
                    previous_index = find_prev_valid(gene_sequence, x)
                    if previous_index is not None:
                        gene_sequence[x].delta_r = gene_sequence[x].start1 - gene_sequence[previous_index].end1
                        gene_sequence[x].delta_q = gene_sequence[x].start2 - gene_sequence[previous_index].end2
                        gene_sequence[x].delta_x = gene_sequence[x].delta_q - gene_sequence[x].delta_r

                        flag = cutoff_check(gene_sequence[x].delta_r, gene_sequence[x].delta_q, cutoff)
                        if flag:
                            gene_sequence[x].delta_r = "--"
                            gene_sequence[x].delta_q = "--"
                            gene_sequence[x].delta_x = "--"

            # Distance calculations for inversions
            else:
                # Distance calculation for standard inverted mum
                next_index = find_next_valid(gene_sequence, x)
                prev_index = find_prev_valid(gene_sequence, x)
                if next_index is not None and prev_index is not None:
                    if gene_sequence[x].inv_last is True:
                        non_inv = find_prev_not_inverted(gene_sequence, x)
                        if non_inv is not None:
                            gene_sequence[x].delta_r = gene_sequence[x].start1 - gene_sequence[prev_index].end1
                            gene_sequence[x].delta_q = gene_sequence[x].start2 - gene_sequence[non_inv].end2
                            gene_sequence[x].delta_x = gene_sequence[x].delta_r - gene_sequence[x].delta_q
                    # Distance calculation for head of the inversion
                    elif gene_sequence[x].inv_head is True:
                        gene_sequence[x].delta_r = gene_sequence[x].start1 - gene_sequence[prev_index].end1
                        gene_sequence[x].delta_q = gene_sequence[x].start2 - gene_sequence[next_index].end2
                        gene_sequence[x].delta_x = gene_sequence[x].delta_q - gene_sequence[x].delta_r
                    else:
                        gene_sequence[x].delta_r = gene_sequence[x].start1 - gene_sequence[prev_index].end1
                        gene_sequence[x].delta_q = gene_sequence[x].start2 - gene_sequence[next_index].end2
                        gene_sequence[x].delta_x = gene_sequence[x].delta_q - gene_sequence[x].delta_r

                    flag = cutoff_check(gene_sequence[x].delta_r, gene_sequence[x].delta_q, cutoff)
                    if flag:
                        gene_sequence[x].delta_r = "--"
                        gene_sequence[x].delta_q = "--"
                        gene_sequence[x].delta_x = "--"


def display_help():
    """
    Function to display help information for the user.
    This is displayed when there are no inputted arguments or when the user calls the
    script with a "-h" flag as an argument.

    :return: None
    """
    help_info = ""
    print("Help Information for the ExpanGe program:")
    print(help_info)


def main(argv):
    """
    The main function for the script. Calls the helper functions to process the inputted data.
    :param argv: The argument list
    :return: None
    """
    file_name = None
    output_name = "ExpanGeOutput.coords"
    sequence = []
    argument_list = argv
    cutoff = 1000000
    options = "hi:o:g:c:"
    long_options = ["help", "input", "output", "geneMap", "cutoff"]
    gene_map = GeneMap()
    gene_map_filename = None

    try:
        flags, values = getopt.getopt(argument_list, options, long_options)

        for current_flag, current_value in flags:
            if current_flag in ("-h", "--help"):
                display_help()
                return None
            elif current_flag in ("-i", "--input"):
                file_name = current_value
            elif current_flag in ("-g", "--geneMap"):
                gene_map_filename = current_value
            elif current_flag in ("-o", "--output"):
                output_name = current_value
            elif current_flag in ("-c", "--cutoff"):
                cutoff = int(current_value)

    except getopt.error as err:
        print(str(err))

    if file_name is None or file_name == "":
        display_help()
        return None

    if gene_map_filename is None or gene_map_filename == "":
        print("No gene map data inputted.")
        print("Gene map must be created and inputted in the command line for proper calculations.")
        print("Help information displayed below.\n")
        display_help()
        return None

    gene_map.create_map(gene_map_filename)
    if gene_map.map is None:
        return None

    # Opens the file and reads through each line.
    # For loop will store the data on each line in a database for later usage
    try:
        data_file = open(file_name, 'r')
    except FileNotFoundError:
        print("Error: The inputted file was not found or could not be opened.")
        print("EXITING PROGRAM")
        return None

    lines = data_file.readlines()
    start_positions = set()
    count = 0

    for line in lines:
        temp_gene = Gene()

        # setting the values for the gene class
        data_list = line.split()
        if len(data_list) != 0:
            temp_gene.start1 = int(data_list[0])
            temp_gene.end1 = int(data_list[1])
            temp_gene.start2 = int(data_list[2])
            temp_gene.end2 = int(data_list[3])
            temp_gene.length1 = int(data_list[4])
            temp_gene.length2 = int(data_list[5])
            temp_gene.IDY = data_list[6]
            temp_gene.ref_chr = data_list[7]
            temp_gene.query_chr = data_list[8]
            temp_gene.inv_count = "--"
            temp_gene.delta_r = "--"
            temp_gene.delta_q = "--"
            temp_gene.delta_x = "--"

        # checking if a transposition has taken place, if there has been, ignore the line
        # To Do: standarized file input, it looks like mums columns are differ based on output

        same_flag = gene_map.same_chromosome(ref_chrom=temp_gene.ref_chr,
                                             ref_pos=temp_gene.start1,
                                             query_chrom=temp_gene.query_chr,
                                             query_pos=temp_gene.start2)
        if same_flag is False:
            temp_gene.ignore = True

        # checking for multiples of genes
        if temp_gene.start1 in start_positions:
            temp_gene.ignore = True
        else:
            start_positions.add(temp_gene.start1)

        sequence.append(temp_gene)

        count = count + 1

    identify_inversions(sequence)
    calculate_distances(sequence, cutoff)

    # Open the output file for writing, create it if it does not exist
    output = open(output_name, "w+")

    # Declare and initialize a pandas dataframe
    fields = ["start1", "end1", "start2", "end2", "length1", "length2", "IDY", "ref_chr", "query_chr", "delta_r", "delta_q", "delta_x", "inv_count"]
    dataframe = pd.DataFrame([vars(f) for f in sequence], columns=fields)
    dataframe.rename(columns={"delta_r": "Delta R", "delta_q": "Delta Q", "delta_x": "Delta X", "inv_count": "Inversion Count"})
    dataframe.to_csv(output, sep="\t", index=False, header=True)


if __name__ == "__main__":
    main(sys.argv[1:])
