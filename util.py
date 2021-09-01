"""
Author: Matthew Harper
File: util.py

This file contains a series of utility functions and classes that are used in the ExpanGe program.
The utilities are imported over to the ExpanGe.py file.
"""

from typing import TypeVar, List, Tuple
import datetime

T = TypeVar("T")
Node = TypeVar("Node")


class Node:
    """
    Implementation of a doubly linked list node.
    Do not modify.
    """
    __slots__ = ["value", "next", "prev"]

    def __init__(self, value: T, next: Node = None, prev: Node = None) -> None:
        """
        Construct a doubly linked list node.

        :param value: value held by the Node.
        :param next: reference to the next Node in the linked list.
        :param prev: reference to the previous Node in the linked list.
        :return: None.
        """
        self.next = next
        self.prev = prev
        self.value = value

    def __repr__(self) -> str:
        """
        Represents the Node as a string.

        :return: string representation of the Node.
        """
        return str(self.value)

    def __str__(self) -> str:
        """
        Represents the Node as a string.

        :return: string representation of the Node.
        """
        return str(self.value)


class DoublyLinkedList:
    __slots__ = ["head", "tail", "size"]

    def __init__(self) -> None:
        """
        Construct an empty doubly linked list.

        :return: None.
        """
        self.head = self.tail = None
        self.size = 0

    def __repr__(self) -> str:
        """
        Represent the DLL as a string.

        :return: string representation of the DLL.
        """
        result = ""
        node = self.head
        while node is not None:
            result += str(node)
            if node.next is not None:
                result += " <-> "
            node = node.next
        return result

    def __str__(self) -> str:
        """
        Represent the DLL as a string.

        :return: string representation of the DLL.
        """
        return repr(self)

    def push(self, val: T, back: bool = True) -> None:
        """
        Create Node containing `val` and add to back (or front) of DLL. Increment size by one.

        :param val: value to be added to the DLL.
        :param back: if True, add Node containing value to back (tail-end) of DLL;
            if False, add to front (head-end).
        :return: None.
        """

        new_node = Node(val)
        self.size = self.size + 1

        if back:
            new_node.prev = self.tail

            if self.tail is None:
                self.head = new_node
                self.tail = new_node
                new_node.next = None
            else:
                self.tail.next = new_node
                new_node.next = None
                self.tail = new_node

        else:
            new_node.next = self.head
            new_node.prev = None

            if self.head is None:
                self.head = new_node
                self.tail = new_node
                new_node.prev = None

            if self.head is not None:
                self.head.prev = new_node
                new_node.prev = None

            self.head = new_node

        pass

    def pop(self, back: bool = True) -> None:
        """
        Remove Node from back (or front) of DLL. Decrement size by 1. If DLL is empty, do nothing.

        :param back: if True, remove Node from (tail-end) of DLL;
            if False, remove from front (head-end).
        :return: None.
        """

        if self.size == 0:
            pass
        elif self.size == 1:
            self.tail = None
            self.head = None
        elif back:
            self.tail = self.tail.prev
            self.tail.next = None
        else:
            self.head = self.head.next
            self.head.prev = None

        self.size = self.size - 1


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
        """
        # The distance between the reference nodes (start1 - end1)
        self.delta_r = None
        # The distance between the query nodes (start2 - end2)
        self.delta_q = None
        # The difference of delta_r and delta_q (delta_r - delta_q)
        self.delta_x = None
        # The number of inversion sections at this point
        self.inv_count = None
        """
        # List of the delta values calculated
        self.delta_list = []
