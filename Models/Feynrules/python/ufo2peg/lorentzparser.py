"""
AST visitor class to extract left / right couplings from Lorentz structures
"""
import ast

def parse_lorentz(expr):
    result = LorentzParser().parse(expr)
    return result


class LorentzParserException(Exception):
    """Base class for all LorentzParser exceptions."""


def cleaner(text):
    """Convert LorentzParser output to a (left,right) tuple of prefactors"""
    if 'm' in text and 'p' in text:
        raise LorentzParserException('Both projections in one term: %s' % text)

    if text in [['g'],['i']]:
        return (1,1)
    else:
        try:
            text.remove('g')
        except ValueError:
            pass
        try:
            text.remove('i')
        except ValueError:
            pass

    def helper(text,char):
        text.remove(char)
        try:
            return int(''.join(text))
        except ValueError:
            text.append('1')
            return int(''.join(text))

    left, right = 0, 0
    if 'm' in text:
        left = helper(text,'m')
    elif 'p' in text:
        right = helper(text,'p')
    elif 'y' in text:
        val = helper(text,'y')
        left = -val
        right = val
    return (left,right)


class LorentzParser(ast.NodeVisitor):
    """Convert UFO Lorentz structures to left / right couplings

    This parser is very sensitive to changes in the way UFO Lorentz
    structures are written out.
    """

    def parse(self,expression):
        from functools import reduce
        self.result = [[]]
        tree = ast.parse(expression)
        #print ('---\n',ast.dump(tree),'\n---')
        cleaned = map(cleaner,self.visit(tree))
        return  reduce(lambda a,b: (a[0]+b[0], a[1]+b[1]), cleaned)

    ##################################

    def newterm(self):
        self.result.append([])

    def add(self,a):
        self.result[-1].append(a)

    ##################################

    def visit_Module(self,node):
        self.generic_visit(node)
        return self.result

    def generic_visit(self,node):
        typename = type(node).__name__
        harmless = ['Module','Expr','BinOp','Mult']
        if typename not in harmless:
            raise LorentzParserException('Missing implementation for %s' % typename)
        super(LorentzParser,self).generic_visit(node)

    def visit_Call(self,node):
        subs = { 
            'Gamma' : 'g',
            'Gamma5' : 'y',
            'Identity' : 'i',
            'ProjM' : 'm',
            'ProjP' : 'p',
        }
        try:
            self.add(subs[node.func.id])
        except KeyError:
            err = "Unknown lorentz component %s" % node.func.id
            raise LorentzParserException(err)

    def visit_Num(self,node):
        self.add(str(node.n))

    def visit_Sub(self,node):
        self.newterm()
        self.add('-')

    def visit_Add(self,node):
        self.newterm()
