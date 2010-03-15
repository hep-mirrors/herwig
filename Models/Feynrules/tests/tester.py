#! /usr/bin/env python

import string


def getTemplate(basename):
    f = open('../%s.template' % basename, 'r')
    templateText = f.read()
    f.close()
    return string.Template( templateText )

def writeFile(filename, text):
    f = open(filename,'w')
    f.write(text)
    f.close()

particles = [ { 'name' : 'bla',
                'pdg' : 666666,
                'mass' : 123.45,
                'spin' : 2,
                'q'   : -3 },
              { 'name' : 'blabar',
                'pdg' : -666666,
                'mass' : 123.45,
                'spin' : 2,
                'q'   : 3 } ]

particleT = string.Template("""create /ThePEG/ParticleData $name
setup $name $pdg $name $mass 0 0 0 $q 0 $spin 1
""")


plist = ""

for p in particles:
    plist += particleT.substitute(p)


subs = [ { 'classname'   : 'blablaZ',
           'BaseSpinTag' : 'FFV' },
 #        { 'classname'   : 'FFP',
 #          'BaseSpinTag' : 'FFV' },
         ]

vertexT = getTemplate( 'Vertex.cc' )


referenceT = string.Template("""
static Reference<FeynRulesModel,ThePEG::Helicity::AbstractFFVVertex> interface$vertex
  ("Vertex/$name", // Vertex/FFZ
   "Reference to the FeynRules Model $vertex",
   &FeynRulesModel::$vertex, false, false, true, false);
""")

getterT = string.Template("""
  tAbstractFFVVertexPtr  $invvertex() const {
    return $vertex;
  }
""")

class ModelSubstitutions:
    def __init__(self):
        self.addVertex = ""
        self.istream = "is"
        self.ostream = "os"
        self.refs    = ""
        self.decls   = ""
        self.getters = ""

    def add(self,name):
        vertexname = '%sVertex' % name
        invname    = 'vertex%s' % name
        self.addVertex += '  addVertex(%s);\n' % vertexname
        self.istream   += ' >> %s' % vertexname
        self.ostream   += ' << %s' % vertexname
        self.refs      += referenceT.substitute(vertex=vertexname,name=name)
        self.decls     += '  AbstractFFVVertexPtr %s;\n' % vertexname
        self.getters   += getterT.substitute(vertex=vertexname,invvertex=invname)

    def __call__(self):
        return { 'addVertex' : self.addVertex,
                 'istream'   : self.istream,
                 'ostream'   : self.ostream,
                 'refs'      : self.refs,
                 'decls'     : self.decls,
                 'getters'   : self.getters }


msubs = ModelSubstitutions()

for sub in subs:
    name = sub['classname']
    vertexname = '%sVertex' % name
    writeFile( 'FR%s.cc' % vertexname, vertexT.substitute(sub) )
    msubs.add(name)

modelT_h  = getTemplate( 'Model.h' )
writeFile( 'Model.h' , modelT_h.substitute(msubs()) )

modelT_cc = getTemplate( 'Model.cc' )
writeFile( 'Model.cc' , modelT_cc.substitute(msubs()) )

modelT_in = getTemplate( 'FR.model' )
writeFile( 'FR.model' , modelT_in.substitute(msubs(),plist=plist) )










#print ccT.substitute(subs)



# make model file from template, switchable whether want our SM or FR's SM

# loop over vertices, make vertex files from template

# produce .in files for SM fundamental particles (if needed) + BSM particles from FR

