// -*- C++ -*-
#ifndef HERWIG_VertexBase_H
#define HERWIG_VertexBase_H
//
// This is the declaration of the VertexBase class.

#include <ThePEG/Config/Complex.h>
#include "ThePEG/Interface/Interfaced.h"
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/PDT/WidthGenerator.h>
#include <ThePEG/Helicity/HelicityDefinitions.h>

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 * 
 *  The VertexBase class is the base class for all helicity amplitude
 *  vertices inside Herwig++. In implements the storage of the particles 
 *  which are allowed to interact at the vertex and some simple functions 
 *  which are often needed by the classes which implement the specific 
 *  vertices.
 *
 */
class VertexBase  : public Interfaced {
  
friend ostream & operator<<(ostream &, const VertexBase &);

public:
  
  /**
   * Default constructor and destructor.
   */
  inline VertexBase();
  inline VertexBase(const VertexBase &);
  virtual ~VertexBase();

  /**
   * Constructors for three point vertices.
   */
  inline VertexBase(int,int,int,vector<int>,vector<int>,vector<int>);

  /**
   * Constructors for four  point vertices.
   */
  inline VertexBase(int,int,int,int,vector<int>,vector<int>,
		    vector<int>,vector<int>);

  /**
   * Constructors for five point vertices.
   */
  inline VertexBase(int,int,int,int,int,vector<int>,vector<int>,
		    vector<int>,vector<int>,vector<int>);

  /**
   * n-point constructor.
   */
  inline VertexBase(int);
  
public:
   
  /**
   * Standard functions for writing and reading from persistent streams.
   */
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
public:

  /**
   * Add particles to the lists.
   */
  
  /**
   * Add item to three point list.
   */
  inline void add(int,int,int);

  /**
   * Add item to four point list.
   */
  inline void add(int,int,int,int);

  /**
   * Add item to five  point list.
   */
  inline void add(int,int,int,int,int);

  /**
   * Number of different particle combinations allowed.
   */
  inline unsigned int size();

  /**
   * Get the number of external particles.
   */
  inline int getNpoint();

  /**
   * Is a particle allowed as an incoming particle?
   */
  inline bool incoming(int);

  /**
   * Is a particle allowed as an outgoing particle?
   */
  inline bool outgoing(int);

  /**
   * Get the list of incoming particles.
   */
  inline vector<PDPtr> getIncoming();

  /**
   * Get the list of outgoing particles.
   */
  inline vector<PDPtr> getOutgoing();

  /**
   * Get the coupling.
   */
  inline const Complex & getNorm();

  /**
   * Function to search the list.
   */
  vector<PDPtr> search(int,int);

  /**
   * Is a given combination allowed.
   */
  bool allowed(int,int,int);
  bool allowed(int,int,int,int);
  bool allowed(int,int,int,int,int);

protected:
  
  /**
   * Standard Interfaced virtual functions.
   */
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void doinitrun();
  inline virtual void dofinish();
  
  /**
   * Change all pointers to Interfaced objects to corresponding clones.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  
  /**
   * Return pointers to all Interfaced objects refered to by this.
   */
  inline virtual IVector getReferences();

protected:

  /**
   * Setup the spins of the particles.
   */
  inline void setSpin(int,int,int);
  inline void setSpin(int,int,int,int);
  inline void setSpin(int,int,int,int,int);

  /**
   * Set up the lists of particles.
   */
  inline void setList(vector<int>,vector<int>,vector<int>);
  inline void setList(vector<int>,vector<int>,vector<int>,vector<int>);
  inline void setList(vector<int>,vector<int>,vector<int>,vector<int>,vector<int>);

  /**
   * Set the list of incoming particles.
   */
  inline void setIncoming();

  /**
   * Set the list of outgoing particles.
   */
  inline void setOutgoing();

  /**
   * Set the number of external particles.
   */
  inline void setNpoint(int);

  /**
   * Set the coupling.
   */
  inline void setNorm(const Complex &);

  /**
   * Calculate the propagator for a diagram.
   */
  inline Complex propagator(int, Energy2,tcPDPtr);

  /**
   * Calculate propagator multiplied by coupling.
   */
  inline Complex normPropagator(int, Energy2,tcPDPtr);
    
private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static AbstractClassDescription<Herwig::Helicity::VertexBase> initVertexBase;
  
  /**
   * Private and non-existent assignment operator.
   */
  VertexBase & operator=(const VertexBase &);
  
private:

  /**
   * Storage of the particles.   
   */

  /**
   * PDG codes.
   */
  vector<int> _iparticlea,_iparticleb,_iparticlec,_iparticled,_iparticlee;

  /**
   * Particle data pointers.
   */
  vector<PDPtr> _particlea,_particleb,_particlec,_particled,_particlee;

  /**
   * Spin.
   */
  vector<int> _ispin;
  int _npoint; unsigned int _nsize;

  /**
   * List of allowed incoming particles.
   */
  vector<PDPtr> _inpart; vector <int> _iinpart;

  /**
   * List of allowed outgoing particles.
   */
  vector<PDPtr> _outpart; vector <int> _ioutpart;

  /**
   * The overall coupling.
   */
  Complex _norm;

};
  
/**
 * Output the information on the vertex.
 */
ostream & operator<<(ostream &, const VertexBase &);

}
}


namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of VertexBase.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::VertexBase,1> {
  typedef Interfaced NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::VertexBase>
  : public ClassTraitsBase<Herwig::Helicity::VertexBase> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/Helicity/VertexBase"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwVertex.so"; }

};

}

#include "VertexBase.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VertexBase.tcc"
#endif

#endif /* HERWIG_VertexBase_H */
