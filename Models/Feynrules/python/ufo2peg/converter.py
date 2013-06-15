"""
AST visitor class to convert Python expressions into C++ as used by ThePEG
"""
import ast


def py2cpp(expr):
    """Convert expr to C++ form. Wraps the converter class."""
    result = PyToCpp().parse(expr)
    return result


class PyToCppException(Exception):
    """Base class for all PyToCpp exceptions."""


class PyToCpp(ast.NodeVisitor):
    """Convert Python math expressions into C++.

    Returns a tuple (expr,syms):
    expr -- C++-compatible expression
    syms -- set of all free variables

    Usage:
    >>> expr = '3+2**a*b'
    >>> PyToCpp().parse(expr)
    ('(3.0+(pow(2.0,a)*b))', set(['a', 'b']))

    Note: 
    The converter is currently not generic, it relies on the
    conventions of Feynrules' UFO format on the one hand and ThePEG's
    C++ types on the other.
    """

    def parse(self,expression):
        """Convert expression to C++ format."""
        self.result = []
        self.symbols = set()
        tree = ast.parse(expression)
        #print ast.dump(tree)
        return self.visit(tree)

    ##################################

    def visit_Module(self,node):
        self.generic_visit(node)
        return ''.join(self.result), self.symbols

    def generic_visit(self,node):
        typename = type(node).__name__
        harmless = ['Module','Expr']
        if typename not in harmless:
            raise PyToCppException('Missing implementation for %s' % typename)
        super(PyToCpp,self).generic_visit(node)

    def visit_UnaryOp(self,node):
        self.result.append('(')
        self.visit(node.op)
        self.visit(node.operand)
        self.result.append(')')

    def visit_BinOp(self,node):
        if type(node.op) == type(ast.Pow()):
            return self.pow_node(node)
    
        self.result.append('(')
        self.visit(node.left)
        self.visit(node.op)
        self.visit(node.right)
        self.result.append(')')

    def pow_node(self,node):
        if is_square(node):
            self.result.append('sqr(')
            self.visit(node.left)
            self.result.append(')')
        else:
            self.result.append('pow(')
            self.visit(node.left)
            self.result.append(',')
            self.visit(node.right)
            self.result.append(')')

    def visit_Call(self,node):
        if is_ii(node): 
            self.result.append('ii')
        else:
            self.visit(node.func)
            self.result.append('(')
            for a in node.args:
                self.visit(a)
                self.result.append(',')
            if self.result[-1] == ',':
                del self.result[-1]
            self.result.append(')')
        
    def visit_Attribute(self,node):
        if node.value.id != 'cmath':
            err = "Don't know how to convert %s module." % node.value.id
            raise PyToCppException(err)
        self.result.append(node.attr)

    def visit_Num(self,node):
        # some zeros are encoded as 0j
        if node.n == 0: text = '0.0'
        else:           text = str(float(node.n))
        self.result.append(text)

    def visit_Name(self,node):
        text = str(node.id)
        if text == 'complex': 
            text = 'Complex'
        elif text == 'complexconjugate': 
            text = 'conj'
        elif text not in []:
            self.symbols.add(text)
        self.result.append(text)

    def visit_Mult(self,node):
        self.result.append('*')

    def visit_Add(self,node):
        self.result.append('+')

    def visit_Sub(self,node):
        self.result.append('-')

    def visit_USub(self,node):
        self.result.append('-')

    def visit_UAdd(self,node):
        self.result.append('+')

    def visit_Div(self,node):
        self.result.append('/')

    def visit_Pow(self,node):
        err = "Shold never get here. BinaryOp catches Pow calls."
        raise PyToCppException(err)

### Helpers

def is_square(node):
    """Check if a Pow object is just a square."""
    try:
        return node.right.n == 2.0
    except:
        return False

def is_ii(node):
    """Check if a Call object is just the imaginary unit."""
    try:
        return ( node.func.id == 'complex' 
                 and node.args[0].n == 0
                 and node.args[1].n == 1 )
    except:
        return False


if __name__ == "__main__":
    import doctest
    doctest.testmod()
