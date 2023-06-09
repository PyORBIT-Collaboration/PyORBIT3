import sys
import math
import posix
import os

from orbit.lattice import AccLattice, AccNode, AccActionsContainer

def printAccNodeStructure(accNode, txt = "", txt_shift = ""):
	txt_shift_local = " "
	txt += txt_shift + txt_shift_local + "==== ENTRANCE AccNode = " + accNode.getName()
	txt += " L=" + str(accNode.getLength())
	txt += os.linesep
	for node in accNode.getChildNodes(AccNode.ENTRANCE):
		txt = printAccNodeStructure(node,txt,txt_shift + txt_shift_local*2)
	txt += txt_shift + txt_shift_local*2 + "==== BODY ENTRANCE n parts =" + str(accNode.getnParts())
	txt += os.linesep
	for ind in range(accNode.getnParts()):
		txt += txt_shift + txt_shift_local*3 + "==== body ind="+str(ind)+" L=" + str(accNode.getLength(ind))
		txt += os.linesep
		txt += txt_shift + txt_shift_local*4 + "==== BEFORE"+os.linesep
		nodes = accNode.getChildNodes(AccNode.BODY,ind,AccNode.BEFORE)
		for node in nodes:
			txt = printAccNodeStructure(node,txt,txt_shift + txt_shift_local*4)
		txt += txt_shift + txt_shift_local*4 + "==== AFTER "+os.linesep
		nodes = accNode.getChildNodes(AccNode.BODY,ind,AccNode.AFTER)
		for node in nodes:
			txt = printAccNodeStructure(node,txt,txt_shift + txt_shift_local*4)
	txt += txt_shift + txt_shift_local*2 + "==== BODY EXIT     n parts =" + str(accNode.getnParts())
	txt += txt_shift + txt_shift_local + "==== EXIT AccNode = " + accNode.getName()
	txt += " L=" + str(accNode.getLength())
	txt += os.linesep
	for node in accNode.getChildNodes(AccNode.EXIT):
		txt = printAccNodeStructure(node,txt,txt_shift + txt_shift_local*2)			
	txt += txt_shift + txt_shift_local + "==== END of AccNode = " + accNode.getName()	
	txt += os.linesep
	return txt

print("============== INITIAL NODE STRUCTURE =============")

accNode = AccNode("test_node")
accNode.setnParts(5)

L = 0.
nParts = accNode.getnParts()
for ind in range(nParts):
	accNode.setLength(ind+1.0,ind)
	L += accNode.getLength(ind)
accNode.setLength(L)
for ind in range(nParts):
	accNode.setLength(ind+1.0,ind)
	L += accNode.getLength(ind)

accNode.addChildNode(AccNode("chEntrn0"),AccNode.ENTRANCE)
accNode.addChildNode(AccNode("chEntrn1"),AccNode.ENTRANCE)
accNode.addChildNode(AccNode("chEntrn2"),AccNode.ENTRANCE)

accNode.addChildNode(AccNode("chExit0"),AccNode.EXIT)
accNode.addChildNode(AccNode("chExit1"),AccNode.EXIT)
accNode.addChildNode(AccNode("chExit2"),AccNode.EXIT)

accNode.initialize()

txt = printAccNodeStructure(accNode)
print(txt)
sys.exit()
print("============== REVERSED NODE STRUCTURE =============")
accNode.reverseOrder()	
 
txt = printAccNodeStructure(accNode)
print(txt)

