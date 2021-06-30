#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 14:26:17 2021

@author: mtd
"""

from IntegrateSets import IntegrateSets

# datasets
SWORD=-1 #this will be replaced by an actual read
SoS=-1

# estimates of discharge timeseries and flow law parameters from the reach-level algorithms
Stage1Estimates={
	"FLPE":-1,
	"Q":-1}

# run-time parameters. these would be read from text file
MOIparams = {
	"fUncGage":0.05 #fractional uncertainty for a gage
}

Stage2Estimates=IntegrateSets(SWORD,SoS,MOIparams,Stage1Estimates)
	
print(Stage2Estimates)


