#! /usr/bin/env python

import os,sys,glob,errno,shutil,time     
from optparse import OptionParser   # should be argparse replace options by args

 

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

  
#parser = argparse.ArgumentParser()
#parser.add_argument('--usrinfile', help='usrinfile')
#parser.add_argument('--infile', help='infile')
#parser.add_argument('--definfile', help='definfile')
#parser.add_argument('--formtempdir', help='formtempdir')
#parser.add_argument('--reduction', help='reduction')
#parser.add_argument('--formopt', help='formopt')
#parser.add_argument('--higgseff', help='higgseff')
#parser.add_argument('--makelink', help='makelink', action="store_true")
#parser.add_argument('--makelinkfrom', help='makelinkfrom')
#parser.add_argument('--makelinkto', help='makelinkto')
#args = parser.parse_args()

  
parser = OptionParser()
parser.add_option("-a",'--usrinfile'   ,dest='usrinfile'    , help='usrinfile')
parser.add_option("-b",'--infile'      ,dest='infile'       , help='infile')
parser.add_option("-c",'--definfile'   ,dest='definfile'    , help='definfile')
parser.add_option("-d",'--formtempdir' ,dest='formtempdir'  , help='formtempdir')
parser.add_option("-e",'--reduction'   ,dest='reduction'    , help='reduction')
parser.add_option("-f",'--formopt'     ,dest='formopt'      , help='formopt')
parser.add_option("-g",'--higgseff'    ,dest='higgseff'     , help='higgseff')
parser.add_option("-i",'--makelink'    ,dest='makelink'     , default=False, help='makelink', action="store_true")
parser.add_option("-j",'--makelinkfrom',dest='makelinkfrom' , help='makelinkfrom')
parser.add_option("-k",'--makelinkto'  ,dest='makelinkto'   , help='makelinkto')
(options, args) = parser.parse_args()



if options.makelink and not os.path.islink(options.makelinkto):
   os.symlink(options.makelinkfrom,options.makelinkto)
   exit()
  



if not os.path.isfile(options.usrinfile):
  shutil.copyfile(options.definfile, options.infile)

if os.path.isfile(options.usrinfile):
  shutil.copyfile(options.usrinfile, options.infile)

replacetext(options.infile,"@MODEL@", "model="+options.higgseff)
replacetext(options.infile,"@REDUCTIONPROGRAMS@", options.reduction)
replacetext(options.infile, "@FORMOPT@", options.formopt)
replacetext(options.infile, "@FORMTEMPDIR@", options.formtempdir)

