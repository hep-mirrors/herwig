#! /usr/bin/env python

import os,sys,glob,errno,shutil,argparse,time
 

#  helper to replace all sourceText in fileName with replaceText
def replacetext(fileName, sourceText, replaceText):
    file = open(fileName, "r") 
    text = file.read() 
    file.close()
    file = open(fileName, "w")
    file.write(text.replace(sourceText, replaceText))
    file.close() 

#  helper to build recursivly path
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise   
      
#  helper to find all files of with name in path
def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)      
	  
	  
parser = argparse.ArgumentParser()
parser.add_argument('--infile', help='infile')
parser.add_argument('--definfile', help='definfile')
parser.add_argument('--formtempdir', help='formtempdir')
parser.add_argument('--reduction', help='reduction')
parser.add_argument('--formopt', help='formopt')
parser.add_argument('--higgseff', help='higgseff')
args = parser.parse_args()

if not os.path.isfile(args.infile):
  shutil.copyfile(args.definfile, args.infile)

replacetext(args.infile,"@MODEL@", "model="+args.higgseff)
replacetext(args.infile,"@REDUCTIONPROGRAMS@", args.reduction)
replacetext(args.infile, "@FORMOPT@", args.formopt)
replacetext(args.infile, "@FORMTEMPDIR@", args.formtempdir)

