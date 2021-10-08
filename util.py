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
        # TEMP
        self.tag = None
        # TEMP
        self.scaffold = None
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

