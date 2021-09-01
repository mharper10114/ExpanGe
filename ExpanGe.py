"""
Author: Matthew Harper
File: ExpanGe.py


"""

import sys
import getopt
from util import Node, DoublyLinkedList, Gene


def identify_inversions(linked_list):
    """
    identify_inversions is a function designed to read through the lines in the .coord file and process all the data for later
    use. It identifies where the reversals occur and will set the values for any values that are to be ignored. Uses a
    recursive helper function
    :param linked_list: The linked_list holding all the lines of the file
    :return: None
    """
    head_node = linked_list.head
    inversion_count = 0

    def identify_inversions_helper(root):
        """
        Recursive helper function for the indentify_inversions function.
        Recusively loops through all nodes in the linked_list and records where inversions exist.
        :param root: The current node
        :return: None
        """
        if root.value.ignore is True:
            identify_inversions_helper(root.next)
        else:
            ignore_flag = False
            end_flag = False
            current_start = root.value.start1
            next_node = root.next

            if next_node is None:
                end_flag = True

            while ignore_flag is False and end_flag is False:
                if next_node.value.ignore is False:
                    ignore_flag = True
                else:
                    next_node = next_node.next
                    if next_node is None:
                        end_flag = True

            if end_flag is False:
                next_start = next_node.value.start1

                if next_start < current_start:
                    next_node.value.reversed = True
                    root.value.reversed = True
                    root.value.inv_count = inversion_count
                    next_node.value.inv_count = inversion_count

                identify_inversions_helper(root.next)

    identify_inversions_helper(head_node)


def find_prev_valid(node):
    """
    Helper function that will take in a node and find the closest, previous node that
    is valid. This is used for calculating the distances between valid nodes.
    :param node: A node containing data
    :return: The previous, valid node
    """
    if node.prev is None:
        return None

    prev_node = node.prev

    if prev_node.value.ignore is True:
        prev_node = find_prev_valid(prev_node)

    return prev_node


def find_next_not_inverted(node) -> Node:
    """
    Search function that will find the next node in the data that is not inverted.
    :param node: Current node
    :return: The not inverted node
    """
    if node is not None:
        if node.value.reversed is True:
            return find_next_not_inverted(node.next)
        else:
            return node
    else:
        return None


def find_next_valid(node) -> Node:
    """
    Search function that finds the next valid node in the linked list
    :param node: Current node
    :return: The next valid node
    """
    if node is not None:
        if node.value.ignore is False:
            return node
        else:
            return find_next_valid(node.next)
    else:
        return None


def calculate_distances(linked_list):
    """
    Function to calculate the distances between each gene in the coords file.
    If a gene's previous gene is invalid, it will find the closest valid gene to calculate the
    distance from.
    Utilizes a recursive helper function.
    :param linked_list: The data from the coords file
    :return: None
    """
    def calculate_distances_helper(root):
        if root is None:
            return None

        if root.value.reversed is False:
            prev_node = find_prev_valid(root)
            if prev_node is not None:
                """
                root.value.delta_r = root.value.start1 - prev_node.value.end1
                root.value.delta_q = root.value.start2 - prev_node.value.end2
                root.value.delta_x = root.value.delta_r - root.value.delta_q
                """

                delta_r = root.value.start1 - prev_node.value.end1
                root.value.delta_list.append(delta_r)
                delta_q = root.value.start2 - prev_node.value.end2
                root.value.delta_list.append(delta_q)
                delta_x = delta_r - delta_q
                root.value.delta_list.append(delta_x)

            else:
                root.value.delta_r = "--"
                root.value.delta_q = "--"
                root.value.delta_x = "--"

        else:
            next_node = find_next_valid(root)
            if next_node is not None:
                """
                root.value.delta_r = next_node.value.start1 - root.value.start1
                root.value.delta_q = next_node.value.start2 - root.value.start2
                root.value.delta_x = root.value.delta_r - root.value.delta_q
                """

                delta_r = next_node.value.start1 - root.value.end1
                root.value.delta_list.append(delta_r)
                delta_q = next_node.value.start2 - root.value.end2
                root.value.delta_list.append(delta_q)
                delta_x = delta_r - delta_q
                root.value.delta_list.append(delta_x)
            else:
                root.value.delta_r = "--"
                root.value.delta_q = "--"
                root.value.delta_x = "--"

        calculate_distances_helper(root.next)

    head_node = linked_list.head
    calculate_distances_helper(head_node)


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
    file_name = ""
    output_name = "ExpanGeOutput.coords"
    sequence = DoublyLinkedList()
    argument_list = argv
    options = "hi:o:"
    long_options = ["help", "input", "output"]

    """
    TO DO: Fix flag arguments to not only display help function
    """
    try:
        flags, values = getopt.getopt(argument_list, options, long_options)

        for current_flag, current_value in flags:
            if current_flag in ("-h", "--help"):
                display_help()
            elif current_flag in ("-i", "--input"):
                file_name = current_value
            elif current_flag in ("-o", "--output"):
                output_name = current_value

    except getopt.error as err:
        print(str(err))

    if file_name == "":
        display_help()
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
    header = ""
    count = 0

    for line in lines:
        if count > 20:
            temp_gene = Gene()

            # setting the values for the gene class
            """
            TO DO: Rewrite string splicing to make universal for all .coords files
            TO DO: Switch parsing over to tab delimited file
            """
            temp_gene.start1 = int(line[:7].strip())
            temp_gene.end1 = int(line[8:16].strip())
            temp_gene.start2 = int(line[20:28].strip())
            temp_gene.end2 = int(line[29:37].strip())
            temp_gene.length1 = int(line[41:49].strip())
            temp_gene.length2 = int(line[50:58].strip())
            temp_gene.IDY = line[62:70].strip()
            temp_gene.tag = line[74:96].strip()
            temp_gene.scaffold = line[97:].strip()

            # checking if a transposition has taken place, if there has been, ignore the line
            # if temp_gene.tag[12:] != temp_gene.scaffold:
                # temp_gene.ignore = True

            # checking for multiples of genes
            if temp_gene.start1 in start_positions:
                temp_gene.ignore = True
            else:
                start_positions.add(temp_gene.start1)

            if temp_gene.ignore is True:
                temp_gene.delta_r = "--"
                temp_gene.delta_x = "--"
                temp_gene.delta_q = "--"
                temp_gene.inv_count = "--"

            sequence.push(temp_gene)

        else:
            header = header + line
            count = count + 1

    identify_inversions(sequence)
    calculate_distances(sequence)

    # Open the output file for writing, create it if it does not exist
    output = open(output_name, "w+")

    # Go through each line in the input file and write the data to the output file
    # Add on the additional data calculated in this program at the end of each line
    curr_node = sequence.head
    count = 0

    for line in lines:
        curr_node = sequence.head

        if count > 5:
            delta_r = str(curr_node.value.delta_r)
            delta_q = str(curr_node.value.delta_q)
            delta_x = str(curr_node.value.delta_x)
            inversion_count = str(curr_node.value.inv_count)

            line.rstrip("\n")

            """
            TO DO: Fix output file formatting
            """

            line = line + "\t" + delta_r + " | " + delta_q + " | " + delta_x + " | " + inversion_count + "\n"
            output.write(line)
            '''
            TO DO: Add columns at end of line for calculated data
            '''
            curr_node = curr_node.next
        else:
            output.write(line)

        count = count + 1


if __name__ == "__main__":
    main(sys.argv[1:])
