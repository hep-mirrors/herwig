// -*- C++ -*-
#ifndef HERWIG_VertexBase_H
#define HERWIG_VertexBase_H
//
// This is the declaration of the <!id>VertexBase<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// The <!id>VertexBase<!!id> class is the base class for all helicity amplitude
// vertices inside Herwig++. In implements the storage of the particles which are
// allowed to interact at the vertex and some simple functions which are often
// needed by the classes which implement the specific vertices.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:.html">.h</a>,
// <a href="http:.html">.h</a>.
// 
#include <ThePEG/Config/Complex.h>
#include "ThePEG/Interface/Interfaced.h"
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/PDT/WidthGenerator.h>
#include <ThePEG/Helicity/HelicityDefinitions.h>

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

class VertexBase  : public Interfaced {
  
friend ostream & operator<<(ostream &, const VertexBase &);

public:
  
  // default constructor and destructor
  inline VertexBase();
  inline VertexBase(const VertexBase &);
  virtual ~VertexBase();

  // constructors for three point vertices
  inline VertexBase(int,int,int,vector<int>,vector<int>,vector<int>);

  // constructors for four  point vertices
  inline VertexBase(int,int,int,int,vector<int>,vector<int>,
		    vector<int>,vector<int>);

  // constructors for five point vertices
  inline VertexBase(int,int,int,int,int,vector<int>,vector<int>,
		    vector<int>,vector<int>,vector<int>);

  // n-point constructor
  inline VertexBase(int);
  
public:
  
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
public:

  // add particles to the lists
  // add item to three point list
  inline void add(int,int,int);
  // add item to four  point list
  inline void add(int,int,int,int);
  // add item to five  point list
  inline void add(int,int,int,int,int);

  // number of different particle combinations allowed
  inline unsigned int size();

  // get the number of external particles
  inline int getNpoint();

  // is a particle allowed as an incoming particle
  inline bool incoming(int);

  // is a particle allowed as an outgoing particle
  inline bool outgoing(int);

  // get the list of incoming particles
  inline vector<PDPtr> getIncoming();

  // get the list of outgoing particles
  inline vector<PDPtr> getOutgoing();

  // get the coupling
  inline const Complex & getNorm();

  // function to search the list
  vector<PDPtr> search(int,int);

  // is a given combination allowed
  bool allowed(int,int,int);
  bool allowed(int,int,int,int);
  bool allowed(int,int,int,int,int);

protected:
  
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void doinitrun();
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.
  
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.
  
  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.

protected:

  // setup the spins of the particles
  inline void setSpin(int,int,int);
  inline void setSpin(int,int,int,int);
  inline void setSpin(int,int,int,int,int);

  // set up the lists of particles
  inline void setList(vector<int>,vector<int>,vector<int>);
  inline void setList(vector<int>,vector<int>,vector<int>,vector<int>);
  inline void setList(vector<int>,vector<int>,vector<int>,vector<int>,vector<int>);

  // set the list of incoming particles
  inline void setIncoming();

  // set the list of outgoing particles
  inline void setOutgoing();

  // set the number of external particles
  inline void setNpoint(int);

  // set the coupling
  inline void setNorm(const Complex &);

  // calculate the propagator for a diagram
  inline Complex propagator(int, Energy2,tcPDPtr);

  // calculate propagator multiplied by coupling
  inline Complex normPropagator(int, Energy2,tcPDPtr);

protected:
  
  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods.
    
private:
  
  static ClassDescription<Herwig::Helicity::VertexBase> initVertexBase;
  // Describe a concrete class with persistent data.
  
  VertexBase & operator=(const VertexBase &);
  // Private and non-existent assignment operator.
  
private:

  // storage of the particles
  // pdg codes
  vector<int> _iparticlea,_iparticleb,_iparticlec,_iparticled,_iparticlee;
  // particle data pointers
  vector<PDPtr> _particlea,_particleb,_particlec,_particled,_particlee;
  // spin
  vector<int> _ispin;
  int _npoint; unsigned int _nsize;
  // list of allowed incoming particles
  vector<PDPtr> _inpart; vector <int> _iinpart;
  // list of allowed outgoing particles
  vector<PDPtr> _outpart; vector <int> _ioutpart;
  // the overall coupling
  Complex _norm;

};
  
// output the information on the vertex
ostream & operator<<(ostream &, const VertexBase &);

}
}

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of VertexBase.
template <>
struct BaseClassTrait<Herwig::Helicity::VertexBase,1> {
  typedef Interfaced NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::Helicity::VertexBase>
  : public ClassTraitsBase<Herwig::Helicity::VertexBase> {
  static string className() { return "/Herwig++/Helicity/VertexBase"; }
  // Return the class name.
  static string library() { return "libHwVertex.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "VertexBase.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VertexBase.tcc"
#endif

#endif /* HERWIG_VertexBase_H */
