import re
import os

class Histogram :
    '''storage of a histogram'''
    x  = []
    y  = []
    dx = []
    data_y   = []
    data_dy  = []
    titleTop=""
    caseTop=""
    titleBottom=""
    caseBottom=""
    titleLeft=""
    caseLeft=""
    titleRight=""
    caseRight=""
    hasData = False
    def write(self,file,new,colour,top,bottom,scaleX,scaleY):
        if(new) :
            file.write("NEW FRAME\n")
            file.write("SET ORDER X Y DX DY\n")
            file.write("SET FONT DUPLEX\n")
            file.write("SET WINDOW X 2 12 Y 3.5 9.0\n")
            if(scaleX) :
                file.write("SET SCALE X LOG\n")
            else :
                file.write("SET SCALE Y LIN\n")
            if(scaleY) :
                file.write("SET SCALE Y LOG\n")
            else :
                file.write("SET SCALE Y LIN\n")
            if(top) :
                file.write(self.titleTop)
                file.write(self.caseTop)
                file.write(self.titleLeft)
                file.write(self.caseLeft)
                file.write(self.titleRight)
                file.write(self.caseRight)
                out =  'SET LIMITS X %s \t %s \t \n' % (self.x[0]-self.dx[0],self.x[len(self.x)-1]+self.dx[len(self.x)-1])
                file.write(out)
            if bottom :
                file.write(self.titleBottom)
                file.write(self.caseBottom)
                file.write("SET AXIS BOTTOM ON\n")
            else :
                file.write("SET AXIS BOTTOM OFF\n")
        for i in range(0,len(self.x)) :
            out =  '%s \t %s \t %s \n' % (self.x[i],self.y[i],self.dx[i])
            file.write(out)
        file.write("HIST " + colour  + "\n")

    def writeData(self,file):
        if( not self.hasData) : return
        for i in range(0,len(self.x)) :
            out =  '%s \t %s \t %s \t %s \n' % (self.x[i],self.data_y[i],self.dx[i],self.data_dy[i])
            file.write(out)
        file.write("PLOT\n")

    def writeDifference(self,other,file,new,colour) :
        if(new) :
            file.write("SET WINDOW X 2 12 Y 2.8 3.5\n")
            file.write("SET SCALE Y LIN\n")
            file.write(self.titleBottom)
            file.write(self.caseBottom)
            out =  'SET LIMITS X %s \t %s \t \n' % (self.x[0]-self.dx[0],self.x[len(self.x)-1]+self.dx[len(self.x)-1])
            file.write(out)
            file.write("SET AXIS BOTTOM ON\n")
        ndegree = 0
        chi = 0
        for i in range(0,len(self.x)) :
            diff = 0.
            if (self.y[i]!=0. and other.y[i]!=0) :
                diff = (self.y[i]-other.y[i])/(self.y[i]+other.y[i])
                ndegree += 1
                chi += abs(diff)
            out =  '%s \t %s \t %s \n' % (self.x[i],diff,self.dx[i])
            file.write(out)
        file.write("JOIN LEVEL=1 " + colour  + "\n")
        return [ndegree,chi]

    def __add__(self, other):
        newHist = Histogram()
        newHist.x          = self.x         
        newHist.dx         = self.dx         
        newHist.titleTop   = self.titleTop   
        newHist.caseTop    = self.caseTop    
        newHist.titleBottom= self.titleBottom
        newHist.caseBottom = self.caseBottom 
        newHist.titleLeft  = self.titleLeft  
        newHist.caseLeft   = self.caseLeft   
        newHist.titleRight = self.titleRight 
        newHist.caseRight  = self.caseRight
        newHist.hasData    = self.hasData
        newHist.data_y = self.data_y
        newHist.data_dy = self.data_dy
        newHist.y=[]
        for j in range(len(self.y)):
            newHist.y.append(self.y[j]+other.y[j])
            
        return newHist  

    def __subtract__(self, other):
        newHist = Histogram()
        newHist.x          = self.x         
        newHist.dx         = self.dx         
        newHist.titleTop   = self.titleTop   
        newHist.caseTop    = self.caseTop    
        newHist.titleBottom= self.titleBottom
        newHist.caseBottom = self.caseBottom 
        newHist.titleLeft  = self.titleLeft  
        newHist.caseLeft   = self.caseLeft   
        newHist.titleRight = self.titleRight 
        newHist.caseRight  = self.caseRight  
        newHist.hasData    = self.hasData
        newHist.data_y = self.data_y
        newHist.data_dy = self.data_dy
        newHist.y=[]
        for j in range(len(self.y)):
            newHist.y.append(self.y[j]-other.y[j])
        return newHist  

    def __mul__(self, other):
        newHist = Histogram()
        newHist.x          = self.x         
        newHist.dx         = self.dx         
        newHist.titleTop   = self.titleTop   
        newHist.caseTop    = self.caseTop    
        newHist.titleBottom= self.titleBottom
        newHist.caseBottom = self.caseBottom 
        newHist.titleLeft  = self.titleLeft  
        newHist.caseLeft   = self.caseLeft   
        newHist.titleRight = self.titleRight 
        newHist.caseRight  = self.caseRight  
        newHist.hasData    = self.hasData
        newHist.data_y = self.data_y
        newHist.data_dy = self.data_dy
        newHist.y=[]
        for j in range(len(self.y)):
            newHist.y.append(self.y[j]*other)
        return newHist  

    def __div__(self, other):
        newHist = Histogram()
        newHist.x          = self.x         
        newHist.dx         = self.dx         
        newHist.titleTop   = self.titleTop   
        newHist.caseTop    = self.caseTop    
        newHist.titleBottom= self.titleBottom
        newHist.caseBottom = self.caseBottom 
        newHist.titleLeft  = self.titleLeft  
        newHist.caseLeft   = self.caseLeft   
        newHist.titleRight = self.titleRight 
        newHist.caseRight  = self.caseRight  
        newHist.hasData    = self.hasData
        newHist.data_y = self.data_y
        newHist.data_dy = self.data_dy
        newHist.y=[]
        for j in range(len(self.y)):
            newHist.y.append(self.y[j]/other)
        return newHist

def readHistogram(f1,data,nplot=1) :
    test = f1.readline()
    hist=[]
    ncount=0
    for i in range(0,nplot) :
        hist.append(Histogram())
        hist[i].x  = []
        hist[i].y  = []
        hist[i].dx = []
        hist[i].data_y  = []
        hist[i].data_dy = []
        hist[i].titleTop=""
        hist[i].caseTop=""
        hist[i].titleBottom=""
        hist[i].caseBottom=""
        hist[i].titleLeft=""
        hist[i].caseLeft=""
        hist[i].titleRight=""
        hist[i].caseRight=""
        hist[i].hasData = data
    while test :
        test = test.lstrip()
        if (test[0] == "N" or test[0:5] == "SET P" or test[0:6] == "SET AX" or
            test[0:5] == "SET W" or test[0:5] == "SET O" or test[0:5] == "SET F" or test[0:5] == "SET L" or test[0:5] == "SET S" ) :
            pass
        elif (test[0:7] == "TITLE B") :
            hist[0].titleBottom = test
            test = f1.readline()
            hist[0].caseBottom = test
        elif (test[0:7] == "TITLE T") :
            hist[0].titleTop = test
            test = f1.readline()
            hist[0].caseTop = test
        elif( test[0:7] == "TITLE L") :
            hist[0].titleLeft = test
            test = f1.readline()
            hist[0].caseLeft = test
        elif (test[0:7] == "TITLE R") :
            hist[0].titleRight = test
            test = f1.readline()
            hist[0].caseRight = test
        elif (test[0:4] == "HIST") :
            ncount += 1
            if(ncount==nplot) :
                break
        else :
            temp = test.split()
            hist[ncount].x .append(float(temp[0]))
            hist[ncount].y .append(float(temp[1]))
            hist[ncount].dx.append(float(temp[2]))
        test = f1.readline()
    if(data) :
        test = f1.readline()
        while test :
            test = test.lstrip()
            if (test[0:4] == "PLOT") :
                break
            elif (test[0:7] == "TITLE B") :
                hist[0].titleBottom = test
                test = f1.readline()
                hist[0].caseBottom = test
            else :
                temp = test.split()
                hist[0].data_y .append(float(temp[1]))
                hist[0].data_dy.append(float(temp[3]))
            test = f1.readline()
        test = f1.readline()
        while test :
            if (test[0:4] == "JOIN") :
                break
            elif (test[0:7] == "TITLE B") :
                hist[0].titleBottom = test
                test = f1.readline()
                hist[0].caseBottom = test
            test = f1.readline()
    if(nplot==1) :
        return hist[0]
    else :
        for i in range(1,nplot) :
            hist[i].data_y      = hist[0].data_y
            hist[i].data_dy     = hist[0].data_dy
            hist[i].titleTop    = hist[0].titleTop
            hist[i].caseTop     = hist[0].caseTop
            hist[i].titleBottom = hist[0].titleBottom
            hist[i].caseBottom  = hist[0].caseBottom
            hist[i].titleLeft   = hist[0].titleLeft
            hist[i].caseLeft    = hist[0].caseLeft
            hist[i].titleRight  = hist[0].titleRight
            hist[i].caseRight   = hist[0].caseRight
            hist[i].hasData     = hist[0].hasData
        return hist

def compareTopdrawFiles(fname1,fname2,fo,nplot,data=False) :
    f1 = open(fname1)
    f2 = open(fname2)
    totalDegree = 0
    totalChi = 0.
    for i in range(0,nplot) :
        h1 = readHistogram(f1,data)
        h2 = readHistogram(f2,data)
        h1.write(fo,True,"BLACK",True,False,False,False)
        h2.write(fo,False,"RED",True,False,False,False)
        if(h1.hasData) : h1.writeData(fo)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    f1.close()
    f2.close()
    return [totalDegree,totalChi]

class LinePlot :
    '''storage of a histogram'''
    x  = []
    y  = []
    data_y   = []
    data_dy  = []
    titleTop=""
    caseTop=""
    titleBottom=""
    caseBottom=""
    titleLeft=""
    caseLeft=""
    titleRight=""
    caseRight=""
    hasData = False
    def write(self,file,new,colour,top,bottom,scaleX,scaleY):
        if(new) :
            file.write("NEW FRAME\n")
            file.write("SET ORDER X Y DY\n")
            file.write("SET FONT DUPLEX\n")
            file.write("SET WINDOW X 2 12 Y 3.5 9.0\n")
            if(scaleX) :
                file.write("SET SCALE X LOG\n")
            else :
                file.write("SET SCALE X LIN\n")
            if(scaleY) :
                file.write("SET SCALE Y LOG\n")
            else :
                file.write("SET SCALE Y LIN\n")
            if(top) :
                file.write(self.titleTop)
                file.write(self.caseTop)
                file.write(self.titleLeft)
                file.write(self.caseLeft)
                file.write(self.titleRight)
                file.write(self.caseRight)
                out =  'SET LIMITS X %s \t %s \t \n' % (self.x[0],self.x[len(self.x)-1])
                file.write(out)
            if bottom :
                file.write(self.titleBottom)
                file.write(self.caseBottom)
                file.write("SET AXIS BOTTOM ON\n")
            else :
                file.write("SET AXIS BOTTOM OFF\n")
        for i in range(0,len(self.x)) :
            out =  '%s \t %s\n' % (self.x[i],self.y[i])
            file.write(out)
        file.write("JOIN LEVEL=1 " + colour  + "\n")

    def writeData(self,file):
        if( not self.hasData) : return
        for i in range(0,len(self.x)) :
            out =  '%s \t %s \t %s \n' % (self.x[i],self.data_y[i],self.data_dy[i])
            file.write(out)
        file.write("PLOT\n")

    def writeDifference(self,other,file,new,colour) :
        if(new) :
            file.write("SET WINDOW X 2 12 Y 2.6 3.5\n")
            file.write("SET SCALE Y LIN\n")
            file.write(self.titleBottom)
            file.write(self.caseBottom)
            out =  'SET LIMITS X %s \t %s\n' % (self.x[0],self.x[len(self.x)-1])
            file.write(out)
            file.write("SET AXIS BOTTOM ON\n")
        ndegree = 0
        chi = 0
        for i in range(0,len(self.x)) :
            diff = 0.
            if (self.y[i]!=0. and other.y[i]!=0) :
                diff = (self.y[i]-other.y[i])/(self.y[i]+other.y[i])
                ndegree += 1
                chi += abs(diff)
            out =  '%s \t %s \n' % (self.x[i],diff)
            file.write(out)
        file.write("JOIN LEVEL=1 " + colour  + "\n")
        return [ndegree,chi]

    def __add__(self, other):
        newLine = LinePlot()
        newLine.x          = self.x          
        newLine.titleTop   = self.titleTop   
        newLine.caseTop    = self.caseTop    
        newLine.titleBottom= self.titleBottom
        newLine.caseBottom = self.caseBottom 
        newLine.titleLeft  = self.titleLeft  
        newLine.caseLeft   = self.caseLeft   
        newLine.titleRight = self.titleRight 
        newLine.caseRight  = self.caseRight
        newLine.hasData    = self.hasData
        newLine.data_y = self.data_y
        newLine.data_dy = self.data_dy
        newLine.y=[]
        for j in range(len(self.y)):
            newLine.y.append(self.y[j]+other.y[j])
            
        return newLine  

    def __subtract__(self, other):
        newLine = LinePlot()
        newLine.x          = self.x          
        newLine.titleTop   = self.titleTop   
        newLine.caseTop    = self.caseTop    
        newLine.titleBottom= self.titleBottom
        newLine.caseBottom = self.caseBottom 
        newLine.titleLeft  = self.titleLeft  
        newLine.caseLeft   = self.caseLeft   
        newLine.titleRight = self.titleRight 
        newLine.caseRight  = self.caseRight  
        newLine.hasData    = self.hasData
        newLine.data_y = self.data_y
        newLine.data_dy = self.data_dy
        newLine.y=[]
        for j in range(len(self.y)):
            newLine.y.append(self.y[j]-other.y[j])
        return newLine  

    def __mul__(self, other):
        newLine = LinePlot()
        newLine.x          = self.x          
        newLine.titleTop   = self.titleTop   
        newLine.caseTop    = self.caseTop    
        newLine.titleBottom= self.titleBottom
        newLine.caseBottom = self.caseBottom 
        newLine.titleLeft  = self.titleLeft  
        newLine.caseLeft   = self.caseLeft   
        newLine.titleRight = self.titleRight 
        newLine.caseRight  = self.caseRight  
        newLine.hasData    = self.hasData
        newLine.data_y = self.data_y
        newLine.data_dy = self.data_dy
        newLine.y=[]
        for j in range(len(self.y)):
            newLine.y.append(self.y[j]*other)
        return newLine  

    def __div__(self, other):
        newLine = LinePlot()
        newLine.x          = self.x          
        newLine.titleTop   = self.titleTop   
        newLine.caseTop    = self.caseTop    
        newLine.titleBottom= self.titleBottom
        newLine.caseBottom = self.caseBottom 
        newLine.titleLeft  = self.titleLeft  
        newLine.caseLeft   = self.caseLeft   
        newLine.titleRight = self.titleRight 
        newLine.caseRight  = self.caseRight  
        newLine.hasData    = self.hasData
        newLine.data_y = self.data_y
        newLine.data_dy = self.data_dy
        newLine.y=[]
        for j in range(len(self.y)):
            newLine.y.append(self.y[j]/other)
        return newLine

def readLine(f1,data=False) :
    test = f1.readline()
    hist = LinePlot()
    hist.x  = []
    hist.y  = []
    hist.data_y  = []
    hist.data_dy = []
    hist.titleTop=""
    hist.caseTop=""
    hist.titleBottom=""
    hist.caseBottom=""
    hist.titleLeft=""
    hist.caseLeft=""
    hist.titleRight=""
    hist.caseRight=""
    hist.hasData = data
    while test :
        test = test.lstrip()
        if (test[0] == "N" or test[0:5] == "SET P" or test[0:6] == "SET AX" or
            test[0:5] == "SET W" or test[0:5] == "SET O" or test[0:5] == "SET F" or test[0:5] == "SET L" or test[0:5] == "SET S" ) :
            pass
        elif (test[0:7] == "TITLE B") :
            hist.titleBottom = test
            test = f1.readline()
            hist.caseBottom = test
        elif (test[0:7] == "TITLE T") :
            hist.titleTop = test
            test = f1.readline()
            hist.caseTop = test
        elif( test[0:7] == "TITLE L") :
            hist.titleLeft = test
            test = f1.readline()
            hist.caseLeft = test
        elif (test[0:7] == "TITLE R") :
            hist.titleRight = test
            test = f1.readline()
            hist.caseRight = test
        elif (test[0:4] == "JOIN") :
            break
        else :
            temp = test.split()
            hist.x .append(float(temp[0]))
            hist.y .append(float(temp[1]))
        test = f1.readline()
    if(data) :
        test = f1.readline()
        while test :
            test = test.lstrip()
            if (test[0:4] == "PLOT") :
                break
            elif (test[0:7] == "TITLE B") :
                hist.titleBottom = test
                test = f1.readline()
                hist.caseBottom = test
            else :
                temp = test.split()
                hist.data_y .append(float(temp[1]))
                hist.data_dy.append(float(temp[2]))
            test = f1.readline()
        test = f1.readline()
        while test :
            if (test[0:4] == "JOIN") :
                break
            elif (test[0:7] == "TITLE B") :
                hist.titleBottom = test
                test = f1.readline()
                hist.caseBottom = test
            test = f1.readline()
        test = f1.readline()
        while test :
            if (test[0:4] == "JOIN") :
                break
            elif (test[0:7] == "TITLE B") :
                hist.titleBottom = test
                test = f1.readline()
                hist.caseBottom = test
            test = f1.readline()
    return hist
