#! /usr/bin/python

import sys
import os

print "\n----------Init Mesh ------------\n"
os.system("Control-Mesh");

print "\n----------Initial Condition ------------\n"
os.system("Control-Initial");

print "\n----------Forward ------------\n"
os.system("Control-Forward Results/forward.00000.bup 1 20");

print "\n----------Terminal Condition ------------\n"
os.system("Control-Terminal Results/forward.00020.bup 20");

print "\n----------Backward ------------\n"
os.system("Control-Backward Results/backward.00020.bup Results/forward 0 19")
