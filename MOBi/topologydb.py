#!/usr/bin/env python3
#####################################
#
# Filename : topologydb.py
#
# Projectname : diSTruct
#
# Author : Oskar Taubert
#
# Creation Date : Thu 21 Jun 2018 10:19:27 AM CEST
#
# Last Modified : Mon 30 Jul 2018 03:02:34 PM CEST
#
#####################################


class TopologyDB:
    """
    Contains topologies for different molecular building blocks.

    Accessed by letters in corresponding alphabet.
    """

    # TODO redo this bit

    alphabets = set()

    def __init__(self, topDict=None, letterMap=None):
        self.currentAlphabet = None
        self.letterMap = letterMap
        self.topDict = topDict
        return

    def switch(self, alphabet):
        return

    def __getitem__(self, letter):
        buildingBlock = self.letterMap[letter]
        return self.topDict[buildingBlock]
        return
    pass
