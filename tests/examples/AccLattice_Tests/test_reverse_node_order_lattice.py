import sys
import math
import posix
import os
import pytest

from orbit.lattice import AccLattice, AccNode, AccActionsContainer


def printAccNodeStructure(accNode, txt="", txt_shift=""):
    txt_shift_local = " "
    txt += txt_shift + txt_shift_local + "==== ENTRANCE AccNode = " + accNode.getName()
    txt += " L=" + str(accNode.getLength())
    txt += os.linesep
    for node in accNode.getChildNodes(AccNode.ENTRANCE):
        txt = printAccNodeStructure(node, txt, txt_shift + txt_shift_local * 2)
    txt += txt_shift + txt_shift_local * 2 + "==== BODY ENTRANCE n parts =" + str(accNode.getnParts())
    txt += os.linesep
    for ind in range(accNode.getnParts()):
        txt += txt_shift + txt_shift_local * 3 + "==== body ind=" + str(ind) + " L=" + str(accNode.getLength(ind))
        txt += os.linesep
        txt += txt_shift + txt_shift_local * 4 + "==== BEFORE" + os.linesep
        nodes = accNode.getChildNodes(AccNode.BODY, ind, AccNode.BEFORE)
        for node in nodes:
            txt = printAccNodeStructure(node, txt, txt_shift + txt_shift_local * 4)
        txt += txt_shift + txt_shift_local * 4 + "==== AFTER " + os.linesep
        nodes = accNode.getChildNodes(AccNode.BODY, ind, AccNode.AFTER)
        for node in nodes:
            txt = printAccNodeStructure(node, txt, txt_shift + txt_shift_local * 4)
    txt += txt_shift + txt_shift_local * 2 + "==== BODY EXIT     n parts =" + str(accNode.getnParts())
    txt += txt_shift + txt_shift_local + "==== EXIT AccNode = " + accNode.getName()
    txt += " L=" + str(accNode.getLength())
    txt += os.linesep
    for node in accNode.getChildNodes(AccNode.EXIT):
        txt = printAccNodeStructure(node, txt, txt_shift + txt_shift_local * 2)
    txt += txt_shift + txt_shift_local + "==== END of AccNode = " + accNode.getName()
    txt += os.linesep
    return txt


print("============== INITIAL NODE STRUCTURE =============")
accNode = AccNode("test_node")
accNode.setnParts(5)

L = 0.0
nParts = accNode.getnParts()
for ind in range(nParts):
    accNode.setLength(ind + 1.0, ind)
    L += accNode.getLength(ind)
accNode.setLength(L)
for ind in range(nParts):
    accNode.setLength(ind + 1.0, ind)
    L += accNode.getLength(ind)

accNode.addChildNode(AccNode("chEntrn0"), AccNode.ENTRANCE)
accNode.addChildNode(AccNode("chEntrn1"), AccNode.ENTRANCE)
accNode.addChildNode(AccNode("chEntrn2"), AccNode.ENTRANCE)

accNode.addChildNode(AccNode("chExit0"), AccNode.EXIT)
accNode.addChildNode(AccNode("chExit1"), AccNode.EXIT)
accNode.addChildNode(AccNode("chExit2"), AccNode.EXIT)

accNode.initialize()

initial_node_structure = printAccNodeStructure(accNode)
print(initial_node_structure)


def test_initial_node_structure():
    expected = """ ==== ENTRANCE AccNode = test_node L=15.0
   ==== ENTRANCE AccNode = chEntrn0 L=0.0
    ==== BODY ENTRANCE n parts =1
     ==== body ind=0 L=0.0
      ==== BEFORE
      ==== AFTER 
    ==== BODY EXIT     n parts =1   ==== EXIT AccNode = chEntrn0 L=0.0
   ==== END of AccNode = chEntrn0
   ==== ENTRANCE AccNode = chEntrn1 L=0.0
    ==== BODY ENTRANCE n parts =1
     ==== body ind=0 L=0.0
      ==== BEFORE
      ==== AFTER 
    ==== BODY EXIT     n parts =1   ==== EXIT AccNode = chEntrn1 L=0.0
   ==== END of AccNode = chEntrn1
   ==== ENTRANCE AccNode = chEntrn2 L=0.0
    ==== BODY ENTRANCE n parts =1
     ==== body ind=0 L=0.0
      ==== BEFORE
      ==== AFTER 
    ==== BODY EXIT     n parts =1   ==== EXIT AccNode = chEntrn2 L=0.0
   ==== END of AccNode = chEntrn2
  ==== BODY ENTRANCE n parts =5
   ==== body ind=0 L=1.0
    ==== BEFORE
    ==== AFTER 
   ==== body ind=1 L=2.0
    ==== BEFORE
    ==== AFTER 
   ==== body ind=2 L=3.0
    ==== BEFORE
    ==== AFTER 
   ==== body ind=3 L=4.0
    ==== BEFORE
    ==== AFTER 
   ==== body ind=4 L=5.0
    ==== BEFORE
    ==== AFTER 
  ==== BODY EXIT     n parts =5 ==== EXIT AccNode = test_node L=15.0
   ==== ENTRANCE AccNode = chExit0 L=0.0
    ==== BODY ENTRANCE n parts =1
     ==== body ind=0 L=0.0
      ==== BEFORE
      ==== AFTER 
    ==== BODY EXIT     n parts =1   ==== EXIT AccNode = chExit0 L=0.0
   ==== END of AccNode = chExit0
   ==== ENTRANCE AccNode = chExit1 L=0.0
    ==== BODY ENTRANCE n parts =1
     ==== body ind=0 L=0.0
      ==== BEFORE
      ==== AFTER 
    ==== BODY EXIT     n parts =1   ==== EXIT AccNode = chExit1 L=0.0
   ==== END of AccNode = chExit1
   ==== ENTRANCE AccNode = chExit2 L=0.0
    ==== BODY ENTRANCE n parts =1
     ==== body ind=0 L=0.0
      ==== BEFORE
      ==== AFTER 
    ==== BODY EXIT     n parts =1   ==== EXIT AccNode = chExit2 L=0.0
   ==== END of AccNode = chExit2
 ==== END of AccNode = test_node
"""
    assert initial_node_structure == expected


print("============== REVERSED NODE STRUCTURE =============")
accNode.reverseOrder()
reversed_node_structure = printAccNodeStructure(accNode)
print(reversed_node_structure)


def test_reverse_node_structure():
    expected2 = """ ==== ENTRANCE AccNode = test_node L=15.0
   ==== ENTRANCE AccNode = chExit2 L=0.0
    ==== BODY ENTRANCE n parts =1
     ==== body ind=0 L=0.0
      ==== BEFORE
      ==== AFTER 
    ==== BODY EXIT     n parts =1   ==== EXIT AccNode = chExit2 L=0.0
   ==== END of AccNode = chExit2
   ==== ENTRANCE AccNode = chExit1 L=0.0
    ==== BODY ENTRANCE n parts =1
     ==== body ind=0 L=0.0
      ==== BEFORE
      ==== AFTER 
    ==== BODY EXIT     n parts =1   ==== EXIT AccNode = chExit1 L=0.0
   ==== END of AccNode = chExit1
   ==== ENTRANCE AccNode = chExit0 L=0.0
    ==== BODY ENTRANCE n parts =1
     ==== body ind=0 L=0.0
      ==== BEFORE
      ==== AFTER 
    ==== BODY EXIT     n parts =1   ==== EXIT AccNode = chExit0 L=0.0
   ==== END of AccNode = chExit0
  ==== BODY ENTRANCE n parts =5
   ==== body ind=0 L=5.0
    ==== BEFORE
    ==== AFTER 
   ==== body ind=1 L=4.0
    ==== BEFORE
    ==== AFTER 
   ==== body ind=2 L=3.0
    ==== BEFORE
    ==== AFTER 
   ==== body ind=3 L=2.0
    ==== BEFORE
    ==== AFTER 
   ==== body ind=4 L=1.0
    ==== BEFORE
    ==== AFTER 
  ==== BODY EXIT     n parts =5 ==== EXIT AccNode = test_node L=15.0
   ==== ENTRANCE AccNode = chEntrn2 L=0.0
    ==== BODY ENTRANCE n parts =1
     ==== body ind=0 L=0.0
      ==== BEFORE
      ==== AFTER 
    ==== BODY EXIT     n parts =1   ==== EXIT AccNode = chEntrn2 L=0.0
   ==== END of AccNode = chEntrn2
   ==== ENTRANCE AccNode = chEntrn1 L=0.0
    ==== BODY ENTRANCE n parts =1
     ==== body ind=0 L=0.0
      ==== BEFORE
      ==== AFTER 
    ==== BODY EXIT     n parts =1   ==== EXIT AccNode = chEntrn1 L=0.0
   ==== END of AccNode = chEntrn1
   ==== ENTRANCE AccNode = chEntrn0 L=0.0
    ==== BODY ENTRANCE n parts =1
     ==== body ind=0 L=0.0
      ==== BEFORE
      ==== AFTER 
    ==== BODY EXIT     n parts =1   ==== EXIT AccNode = chEntrn0 L=0.0
   ==== END of AccNode = chEntrn0
 ==== END of AccNode = test_node
"""
    assert reversed_node_structure == expected2
