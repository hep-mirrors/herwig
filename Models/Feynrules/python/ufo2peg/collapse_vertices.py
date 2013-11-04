
class CollapsedVertex(object):
    def __init__(self,v):
        self.name = v.name[2:]
        self.particles_list = [v.particles]
        self.color = tuple(v.color)
        self.lorentz = tuple(v.lorentz)
        self.couplings = tuple(v.couplings.iteritems())

    def __hash__(self):
        return hash((self.color, self.lorentz, self.couplings))

    def __eq__(self,other):
        return ( (self.color, self.lorentz, self.couplings) 
                 == (other.color, other.lorentz, other.couplings) )



def collapse_vertices(vs):
    """Collect particle lists with common structure from the Feynrules vertices."""
    vertices = {}
    for v in vs:
        cv = CollapsedVertex(v)
        if cv in vertices:
            vertices[cv].particles_list += cv.particles_list
            vertices[cv].name += '_%s' % cv.name
        else:
            vertices[cv] = cv
    return vertices.keys()
