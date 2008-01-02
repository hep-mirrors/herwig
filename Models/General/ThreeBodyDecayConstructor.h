// -*- C++ -*-
#ifndef HERWIG_ThreeBodyDecayConstructor_H
#define HERWIG_ThreeBodyDecayConstructor_H
//
// This is the declaration of the ThreeBodyDecayConstructor class.
//

#include "NBodyDecayConstructorBase.h"
#include "ThePEG/Helicity/Vertex/VertexBase.h"
#include "TBDiagram.h"
#include "Herwig++/Decay/General/GeneralThreeBodyDecayer.fh"
#include "ThreeBodyDecayConstructor.fh"

namespace Herwig {
using namespace ThePEG;
using Helicity::VertexType;
using Helicity::VertexBasePtr;

/**
 * A two body decay mode which is a prototype for the 
 * three body mode
 */
struct TwoBodyPrototype {

  /**
   *  Constructor
   */
  TwoBodyPrototype(tPDPtr in, tPDPair out, VertexBasePtr v) :
    incoming(in), outgoing(out), vertex(v) {}

  /**
   *  Incoming particle
   */
  tPDPtr incoming;

  /**
   *  Outgoing particles
   */
  tPDPair outgoing;

  /**
   *  The vertex for the interaction
   */
  VertexBasePtr vertex;
};

/**
 *  A struct to order the particles in the same way as in the DecayMode's
 */
struct ParticleOrdering {
  bool operator()(PDPtr p1, PDPtr p2) {
    return abs(p1->id()) > abs(p2->id()) ||
      ( abs(p1->id()) == abs(p2->id()) && p1->id() > p2->id() ) ||
      ( p1->id() == p2->id() && p1->fullName() > p2->fullName() );
  }
};

/**
 * A set of ParticleData objects ordered as for the DecayMode's
 */
typedef multiset<PDPtr,ParticleOrdering> OrderedParticles;

/**
 * The ThreeBodyDecayConstructor class inherits from the dummy base class
 * NBodyDecayConstructorBase and implements the necessary functions in
 * order to create the 3 body decaymodes for a given set of vertices
 * stored in a Model class.
 *
 * @see \ref ThreeBodyDecayConstructorInterfaces "The interfaces"
 * defined for ThreeBodyDecayConstructor.
 * @see NBodyDecayConstructor
 */
class ThreeBodyDecayConstructor: public NBodyDecayConstructorBase {

public:

  /**
   * The default constructor.
   */
  inline ThreeBodyDecayConstructor();

  /**
   * Function used to determine allowed decaymodes, to be implemented
   * in derived class.
   *@param part vector of ParticleData pointers containing particles in model
   */
  virtual void DecayList(const vector<PDPtr> & part);

protected:


  /**
   * Create the two body prototypes for the decays
   * @param inpart Incoming particle 
   * @param vert The vertex to create decays for
   * @param ilist Which list to search
   * @return A vector a decay modes
   */
  vector<TwoBodyPrototype> createPrototypes(tPDPtr inpart, VertexBasePtr vert,
					unsigned int ilist);

  /**
   * Expand the two body prototype to get the possible
   * threebody diagrams
   * @param proto The two body prototype
   * @param vert The vertex to create decays for
   * @param ilist Which list to search
   */
  vector<TBDiagram> expandPrototype(TwoBodyPrototype proto, VertexBasePtr vert,
				    unsigned int ilist);

  /**
   *  Create the decayer
   * @param The diagrams for the decay
   */
  GeneralThreeBodyDecayerPtr createDecayer(const vector<TBDiagram> & diagrams) const;

  /**
   * Contruct the classname and object name for the Decayer
   * @param incoming The incoming particle
   * @param outgoing The decay products
   * @param objname a string containing the default path of the Decayer object
   */  
  string DecayerClassName(tcPDPtr incoming, const OrderedParticles & outgoing, 
			  string & objname) const;

  /**
   *  Create the DecayMode from the diagrams
   * @param diagrams The diagrams
   */
  void createDecayMode(const vector<TBDiagram> & diagrams);

  /**
   * Get the correct colour factor matrix.
   * @param extpart Vector of external ParticleData pointers
   * @param ncf Set the number of colourflows.
   */
  vector<DVector> getColourFactors(tcPDPtr incoming, const OrderedParticles & outgoing, 
				   unsigned int & ncf) const;

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ThreeBodyDecayConstructor> initThreeBodyDecayConstructor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ThreeBodyDecayConstructor & operator=(const ThreeBodyDecayConstructor &);

private:

  /**
   *  Whether or not to remove on-shell diagrams
   */
  bool _removeOnShell;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ThreeBodyDecayConstructor. */
template <>
struct BaseClassTrait<Herwig::ThreeBodyDecayConstructor,1> {
  /** Typedef of the first base class of ThreeBodyDecayConstructor. */
  typedef Herwig::NBodyDecayConstructorBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ThreeBodyDecayConstructor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ThreeBodyDecayConstructor>
  : public ClassTraitsBase<Herwig::ThreeBodyDecayConstructor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ThreeBodyDecayConstructor"; }
};

/** @endcond */

}

#include "ThreeBodyDecayConstructor.icc"

#endif /* HERWIG_ThreeBodyDecayConstructor_H */
