#! /usr/bin/env python

import os, os.path, sys
import glob


def run():
    xml_files = glob.glob("*integrationBins*-grids.xml")
    text=""
    for file in xml_files:
      xml = open(file, "r")
      text+=xml.read() 
      xml.close()  
    print text.replace("</Grids>\n","").replace("<Grids>\n",""),"</Grids>"
print "<Grids>"
run() 
