"""
Author: Matthew Harper
File: util.py

This file contains a series of utility functions and classes that are used in the ExpanGe program.
The utilities are imported over to the ExpanGe.py file.
"""

from typing import TypeVar, List, Tuple
import datetime


class Gene:
    """
    Gene class that represents each line in the inputted coords file
    It contains the start and end points for both genes as well as the length
    and the scaffold it is contained on.
    """
    def __init__(self):

        # The data from each line in the coords file

        # The start position of the reference gene
        self.start1 = None
        # The end position of the reference gene
        self.end1 = None
        # The start position of the query gene
        self.start2 = None
        # The end position of the query gene
        self.end2 = None
        # The length of the reference gene
        self.length1 = None
        # The length of the query gene
        self.length2 = None
        # TEMP
        self.IDY = None
        # reference chromosome
        self.ref_chr = None
        # Query chromosome
        self.query_chr = None
        # Flag to determine if the gene is inverted
        self.reversed = False
        # Flag to determine whether the software will ignore this gene
        self.ignore = False
        # The distance between the reference nodes (start1 - end1)
        self.delta_r = None
        # The distance between the query nodes (start2 - end2)
        self.delta_q = None
        # The difference of delta_r and delta_q (delta_r - delta_q)
        self.delta_x = None
        # The number of inversion sections at this point
        self.inv_count = None
        # Inversion head flag
        self.inv_head = False
        # Inversion tail flag
        self.inv_tail = False
        # Last inversion flag
        self.inv_last = False


class GeneMap:
    """
    GeneMap class code
    A utility class used to create a gene map that identifies which chromosome a mum in the
    sequence exists on. If the reference mum and query mum exist on different chromosomes
    then the line is marked as a transposition and ignored in calculations

    Takes as input file in the following format:

    Example:
    #ref_chrom  ref_len  query_chrom query_len
    chr1    2000    ChR1A   3000

    More:
    ref_chrom :: str name given by user can be anything but must exactly match the name that appears in mum file
    ref_len :: int length of sequence/contig/chromosome
    query_chrom :: str name given by user can be anything but must exactly match the name that appears in mum file
    query_len :: int  length of sequence/contig/chromsome



    """
    def __init__(self):
        # Stores a dictionary with the chromosome number being the key and
        # the value being a tuple with the end positions for the chromosomes
        # on the reference and the query
        self.map = None

    def create_map(self, filename):
        """
        Function to create and the store the dictionary of the chromosome sizes
        :param filename: The filename of the gene map data
        :return: None
        """
        try:
            map_file = open(filename, 'r')
        except FileNotFoundError:
            print("Error: The inputted file was not found or could not be opened.")
            print("EXITING PROGRAM")
            return None

        lines = map_file.readlines()
        temp_map = {}
        count = 0 # looks like this being done skip a header. We can do this by slicing lines  as well. But I don't want
                  # to mess with this too much just yet.

        for line in lines:
            if count != 0:
                #chrom = "Chr" + str(count)
                #chrom will go unused  in this version
                line_data = line.split()
                temp_ref_chrom = line_data[0]
                temp_ref = int(line_data[1])
                temp_query_chrom = line_data[2]
                temp_query = int(line_data[3])

                tup = (temp_ref, temp_query_chrom, temp_query)

                temp_map[temp_ref_chrom] = tup

            count = count + 1

        self.map = temp_map

    def same_chromosome(self, ref_chrom, ref_pos, query_chrom, query_pos):
        """
        Comparison function to determine if the reference mum and the query mum exist on
        the same chromosome number. Used to determine if transposition has occurred.
        :param ref_pos: The position of the reference mum
        :param query_pos: The position of the query mum
        :return: True or False
        """
        return_bool = False
        if ref_chrom in self.map:
            ref_len, expected_query_chrom, query_len = self.map[ref_chrom]
            if query_chrom == expected_query_chrom:
                if ref_pos < ref_len:
                    if query_pos < query_len:
                        return_bool = True
        return return_bool




