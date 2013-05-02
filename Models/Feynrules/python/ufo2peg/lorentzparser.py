"""
AST visitor class to extract left / right couplings from Lorentz structures
"""
import ast


def parse_lorentz(expr):
    result = LorentzParser().parse(expr)
    return result


class LorentzParserException(Exception):
    """Base class for all LorentzParser exceptions."""


class LorentzParser(ast.NodeVisitor):
    """Convert UFO Lorentz structures to left / right couplings"""

    def parse(self,expression):
        self.result = []
        self.symbols = set()
        tree = ast.parse(expression)
        print ast.dump(tree)
#        return self.visit(tree)

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


if __name__ == "__main__":
    import doctest
    doctest.testmod()
