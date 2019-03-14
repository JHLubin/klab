#!/usr/bin/python
from pyrosetta import PyMOLMover

def pm():
	""" PyMOLMover('127.0.0.1', port=65002) """ 
	return PyMOLMover('127.0.0.1', 65002)