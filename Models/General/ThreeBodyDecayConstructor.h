// -*- C++ -*-
#ifndef HERWIG_ThreeBodyDecayConstructor_H
#define HERWIG_ThreeBodyDecayConstructor_H
//
// This is the declaration of the ThreeBodyDecayConstructor class.
//

#include "NBodyDecayConstructorBase.h"
#include "ThePEG/Helicity/Vertex/VertexBase.h"
#include "TBDiagram.h"
#include "PrototypeVertex.h"
#include "Herwig/Decay/General/GeneralThreeBodyDecayer.fh"

namespace Herwig {
using namespace ThePEG;

using Helicity::VertexBasePtr;

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
    includeIntermediatePhotons_(false),
    interOpt_(0), widthOpt_(1), weakMassCut_(-GeV),
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
   * Create the decayer
   * @param diagrams The diagrams for the decay
   * @param inter Option for intermediates
   */
  GeneralThreeBodyDecayerPtr createDecayer(vector<TBDiagram> & diagrams, 
					   bool inter,double symfac) const;

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
  virtual void createDecayMode(vector<NBDiagram> & mode,
			       bool possibleOnShell,
			       double symfac);

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
   *  Option for intermediate photons
   */
  bool includeIntermediatePhotons_;

  /**
   *  Option for the inclusion of intermediates
   */
  unsigned int interOpt_;

  /**
   *  How to treat the widths of the intermediate particles
   */
  unsigned int widthOpt_;

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
