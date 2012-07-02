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

namespace Herwig {
using namespace ThePEG;

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
  ThreeBodyDecayConstructor() : 
    _removeOnShell(1), _includeTopOnShell(false), _interopt(0), _widthopt(1), 
    _minReleaseFraction(1e-3), _maxBoson(1), _maxList(1), weakMassCut_(-GeV),
    intOpt_(1), relErr_(1e-2) {}

  /**
   * Function used to determine allowed decaymodes, to be implemented
   * in derived class.
   *@param part vector of ParticleData pointers containing particles in model
   */
  virtual void DecayList(const set<PDPtr> & part);

  /**
   * Number of outgoing lines. Required for correct ordering.
   */
  virtual unsigned int numBodies() const { return 3; }

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
   * Create the decayer
   * @param diagrams The diagrams for the decay
   * @param inter Option for intermediates
   */
  GeneralThreeBodyDecayerPtr createDecayer(vector<TBDiagram> & diagrams, 
					   bool inter) const;

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
   * @param inter Option for intermediates
   */
  void createDecayMode(vector<TBDiagram> & diagrams, bool inter);

  /**
   * Get the correct colour factor matrix.
   * @param incoming The incoming particle
   * @param outgoing The outgoing particles
   * @param diagrams The diagrams
   * @param ncf Set the number of colourflows.
   */
  pair<vector<DVector>,vector<DVector> >
  getColourFactors(tcPDPtr incoming, const OrderedParticles & outgoing, 
		   const vector<TBDiagram> & diagrams,
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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
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
  unsigned int _removeOnShell;

  /**
   *  Include on-shell for \f$t\to b W\f$
   */
  bool _includeTopOnShell;

  /**
   *  Option for the inclusion of intermediates
   */
  unsigned int _interopt;

  /**
   *  How to treat the widths of the intermediate particles
   */
  unsigned int _widthopt;

  /**
   * The minimum energy release for a three-body decay as a 
   * fraction of the parent mass
   */
  double _minReleaseFraction;

  /**
   *  Maximum number of EW gauge bosons
   */
  unsigned int _maxBoson;

  /**
   *  Maximum number of particles from the decaying particle list
   */
  unsigned int _maxList;

  /**
   *  Excluded Vertices
   */
  vector<VertexBasePtr> excludedVector_;

  /**
   *  Excluded Vertices
   */
  set<VertexBasePtr> excludedSet_;

  /**
   *  Cut off or decays via the weak current
   */
  Energy weakMassCut_;

  /**
   *  Option for the integration to get the partial width
   */
  unsigned int intOpt_;

  /**
   *  Relative error for partial width integration
   */
  double relErr_;
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

#endif /* HERWIG_ThreeBodyDecayConstructor_H */
