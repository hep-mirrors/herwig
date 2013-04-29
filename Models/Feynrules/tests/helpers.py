from string import Template
"""
Helper functions for the Herwig++ Feynrules converter
"""


class CheckUnique:
    """Uniqueness checker.
    
    An object of this class remembers the value it was called with first.
    Any subsequent call to it will only succeed if the same value is passed again.
    For example,

    >>> f = CheckUnique()
    >>> f(5)
    >>> f(5)
    >>> f(4)
    Traceback (most recent call last):
        ...
    AssertionError
    """
    def __init__(self):
        self.val = None

    def __call__(self,val):
        """Store value on first call, then compare."""
        if self.val is None:
            self.val = val
        else:
            assert( val == self.val )



def is_number(s):
    """Check if a value is a number."""
    try:
        float(s)
    except ValueError:
        return False
    return True



def getTemplate(name):
    """Create a template from a file."""
    with open('../{}.template'.format(name), 'r') as f:
        templateText = f.read()
    return Template( templateText )


def writeFile(filename, text):
    """Write text to a filename."""
    with open(filename,'w') as f:
        f.write(text)



def get_lorentztag(spin):
    """Produce a ThePEG spin tag for the given numeric FR spins."""
    spins = { 1 : 'S', 2 : 'F', 3 : 'V', -1 : 'U', 5 : 'T' }
    result = [ spins[s] for s in spin ]

    def spinsort(a,b):
        """Helper function for ThePEG's FVST spin tag ordering."""
        if a == b: return 0
        for letter in 'FVST':
            if a == letter: return -1
            if b == letter: return  1

    result = sorted(result, cmp=spinsort)
    return ''.join(result)



def banner():
    return """\
===============================================================================================================
______                  ______        _                 __ _   _                        _                      
|  ___|                 | ___ \      | |               / /| | | |                      (_)          _      _   
| |_  ___  _   _  _ __  | |_/ /_   _ | |  ___  ___    / / | |_| |  ___  _ __ __      __ _   __ _  _| |_  _| |_ 
|  _|/ _ \| | | || \_ \ |    /| | | || | / _ \/ __|  / /  |  _  | / _ \| \__|\ \ /\ / /| | / _` ||_   _||_   _|
| | |  __/| |_| || | | || |\ \| |_| || ||  __/\__ \ / /   | | | ||  __/| |    \ V  V / | || (_| |  |_|    |_|  
\_|  \___| \__, ||_| |_|\_| \_|\__,_||_| \___||___//_/    \_| |_/ \___||_|     \_/\_/  |_| \__, |              
            __/ |                                                                           __/ |              
           |___/                                                                           |___/               
===============================================================================================================
generating model/vertex/.model/.in files
please be patient!
===============================================================================================================
"""




#################### ??? #######################


# function that replaces ** with pow(,): used in PyMathToThePEGMath below
def powstring(stringrep,power):
    #if(stringrep[0] == '('):
    #    return 'pow(' + stringrep[0:len(stringrep)-2] + stringrep[len(stringrep)-3:len(stringrep)].replace('**','') + ',' + power + ')'
    #else:
    return 'pow(' + stringrep[0:len(stringrep)-3] + stringrep[len(stringrep)-3:len(stringrep)+1].replace('**','') + ',' + power + ')'

def BrackParams(stringin):
    # to take care of the ** powers over brackets (), append new "parameters" of type (xxxx)
    #	th = acos(1/sqrt(1 + (-pow(MM1,2) + pow(MM2,2) + sqrt(4*pow(MM12,4) + (pow(MM1,2) - pow(MM2,2))**2))**2/(4.*pow(MM12,4))));
    # start by finding the positions of the left and right brackets 
    stringbrack = stringin
    d_string_length = len(stringbrack)
    brack_pos_left = []
    brack_pos_right = []
    for ll in range(0,d_string_length):
        if(stringbrack[ll] is '('):
            brack_pos_left.append(ll)
        if(stringbrack[ll] is ')'):
            brack_pos_right.append(ll)

    paramsin_brack = []
    # loop over the left bracket positions, moving left, and count the number of right brackets
    for le in brack_pos_left:
        LRNUM = -1
        # print 'brack_pos_left', le
        # loop over the input string starting from the position of the bracket
        # count the left and right brackets to the right of the left bracket
        # left brackets are negative, right brackets are positive
        # the matched bracket is found when the left + right = LRNUM = 0
        for lm in range(le+1,d_string_length):
            if(lm in brack_pos_left):
                LRNUM = LRNUM - 1
            if(lm in brack_pos_right):
                LRNUM = LRNUM + 1
            if(LRNUM is 0):
                #print 'BRACKET', le, 'matched with', lm, stringbrack
                #print 'appending', stringbrack[le:lm+1]
                paramsin_brack.append(stringbrack[le:lm+1])
                break
            
    # sort the paramsin_brack according to length, from longest to shortest
    # insert to start of string
    paramsin_brack.sort(key=len, reverse=True)
    # for iii in range(0,len(paramsin_brack)):
    #paramsin.insert(0,paramsin_brack[iii])
    return paramsin_brack

# function that converts the vertex expressions in Python math format:
# called twice: once for the parameters, and once for terms in parentheses (...)
def PyMathToThePEGMath(stringin, paramsin):
    stringout = PyMathToThePEGMath_nb(stringin, paramsin)
    paramsbrack = BrackParams(stringout)
    #print 'paramsbrack', paramsbrack
    if(paramsbrack is not []):
        stringreturn =  PyMathToThePEGMath_nb(stringout, paramsbrack)
    return stringreturn

# function that converts the vertex expressions in Python math format to
# format that can be calculated using ThePEG 
def PyMathToThePEGMath_nb(stringin, paramsin):
    
    # define an array that contains the numbers 0-9 in string form
    numbersarray = ['0','1','2','3','4','5','6','7','8','9']

    # define counters and variables used to detect the positions of powers
    ii = 0
    pos = ''
    posnew = ''
    powpos = ''
    pow_ddg = 0
    power = ''
  

    # add '**' to the end of paramsin[ss], the array of given parameters of the model
    for ss in range(0,len(paramsin)):
        paramsin[ss] = paramsin[ss] + '**'

        #print 'in progress', stringin

    #print 'paramsin', paramsin
    powerchange = []
    # loop over the array of the model parameters and search for them in the given mathematical expression
    # each time a new position with the ** notation is found, replace with the C++ pow(,) notation
    for xx in range(0,len(paramsin)):
        # powerchange contains information on the positions
        # of necessary changes. OBSOLETE: for testing purposes only
        powerchange.append([])
        # reset counter for next variable and position variables
        ii = 0
        pos = 0
        posnew = 0
        # save the length of the string at the beginning of the loop for a parameter
        initial_string_length = len(stringin)
        # scan the string from right to left
        while (initial_string_length-ii >= 0):
            # set the new position of the found 
            posnew = stringin.find(paramsin[xx],initial_string_length-ii)
            #print 'param found', posnew, paramsin[xx]
            # if the position is new, do stuff
            if(posnew is not pos and posnew is not -1):
                # get the position of the power
                powpos = posnew + len(paramsin[xx]) 
                # check if power is single or double digit
                # i.e. -> assuming there are no powers beyond "99"
                power = stringin[powpos]
                if(powpos+1 < initial_string_length):
                    if(stringin[powpos+1] in numbersarray):
                        #print 'power is double digit ', (stringin[powpos+1])
                        pow_ddg = 1
                        power = stringin[powpos] + stringin[powpos+1]
                powerchange[xx].append([ posnew, power ])
                # do the replacement of the ** to pow(,)
                stringin = stringin[:posnew] + stringin[posnew:posnew+len(paramsin[xx])+len(power)].replace(paramsin[xx]+power,powstring(paramsin[xx],power)) + stringin[posnew+len(paramsin[xx])+len(power):]
                #print 'in progress', stringin
            # reset position variable for next point in string
            # increment the counter for the position in string
            pos = posnew
            ii += 1
    # do replacements of 'complex'
    # integers multiplying stuff (add the "."
    stringin = stringin.replace('0j','Complex(0,0)')
    stringin = stringin.replace('complex(0,1)','Complex(0,1.)')
    stringin = stringin.replace('complex','Complex')
    stringin = stringin.replace('cmath.pi', 'M_PI')
    stringin = stringin.replace('cmath.', '')
    stringin = stringin.replace('-(','(-1.)*(')
    for nn in range(0,len(numbersarray)):
        numbersnn = numbersarray[nn]
        stringin = stringin.replace(numbersnn +' *', numbersnn +'. *')
        stringin = stringin.replace(numbersnn+'*Complex',numbersnn+'.*Complex')
        
   
    # print 'final string:'
    # print 'final string', stringin, paramsin
    # reset the parameters with ** for next run of function and return
    for ss in range(0,len(paramsin)):
        paramsin[ss] = paramsin[ss].replace('**','')
    return stringin

# function that replaces alphaS (aS)-dependent variables
# with their explicit form which also contains strongCoupling
def aStoStrongCoup(stringin, paramstoreplace, paramstoreplace_expressions):
    #print stringin
    for xx in range(0,len(paramstoreplace)):
        #print paramstoreplace[xx], paramstoreplace_expressions[xx]
        stringout = stringin.replace(paramstoreplace[xx], '(' +  PyMathToThePEGMath(paramstoreplace_expressions[xx],allparams) + ')')
    stringout = stringout.replace('aS', '(sqr(strongCoupling(q2))/(4.0*Constants::pi))')
    #print 'resulting string', stringout
    return stringout


# function that replaces alphaEW (aEW)-dependent variables
# with their explicit form which also contains weakCoupling
def aEWtoWeakCoup(stringin, paramstoreplace, paramstoreplace_expressions):
    #print stringin
    for xx in range(0,len(paramstoreplace)):
        #print paramstoreplace[xx], paramstoreplace_expressions[xx]
        stringout = stringin.replace(paramstoreplace[xx], '(' +  PyMathToThePEGMath(paramstoreplace_expressions[xx],allparams) + ')')
    stringout = stringout.replace('aEWM1', '(1/(sqr(electroMagneticCoupling(q2))/(4.0*Constants::pi)))')
    #print 'resulting string', stringout
    return stringout




if __name__ == "__main__":
    import doctest
    doctest.testmod()
