// -*- C++ -*-
#ifndef THEPEG_FourBodyDecayConstructor_H
#define THEPEG_FourBodyDecayConstructor_H
//
// This is the declaration of the FourBodyDecayConstructor class.
//

#include "NBodyDecayConstructorBase.h"
#include "ThePEG/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Decay/General/GeneralFourBodyDecayer.fh"
#include "TwoBodyPrototype.h"

namespace Herwig {

using namespace ThePEG;
using Helicity::VertexBasePtr;

/**
 * Here is the documentation of the FourBodyDecayConstructor class.
 *
 * @see \ref FourBodyDecayConstructorInterfaces "The interfaces"
 * defined for FourBodyDecayConstructor.
 */
class FourBodyDecayConstructor: public NBodyDecayConstructorBase {

public:

  /**
   * The default constructor.
   */
  FourBodyDecayConstructor() :
    removeOnShell_(1), interopt_(0), widthopt_(1), 
    minReleaseFraction_(1e-3), maxBoson_(0), maxList_(0) {}

  /**
   *  Destructor
   */
  ~FourBodyDecayConstructor();

  /**
   * Function used to determine allowed decaymodes, to be implemented
   * in derived class.
   * @param particles vector of ParticleData pointers containing 
   * particles in model
   */
  virtual void DecayList(const set<PDPtr> & particles);

  /**
   * Number of outgoing lines. Required for correct ordering.
   */
  virtual unsigned int numBodies() const {return 4;}

  /**
   *  Create a decay mode
   */
  void createDecayMode(vector<PrototypeVertexPtr> & diagrams,
		       bool inter);

  /**
   * Create the decayer
   * @param diagrams The diagrams for the decay
   * @param inter Option for intermediates
   */
  GeneralFourBodyDecayerPtr createDecayer(vector<PrototypeVertexPtr> & diagrams,
					  bool inter) const;

  /**
   * Contruct the classname and object name for the Decayer
   * @param incoming The incoming particle
   * @param outgoing The decay products
   * @param objname a string containing the default path of the Decayer object
   */  
  string DecayerClassName(tcPDPtr incoming, const OrderedParticles & outgoing, 
			  string & objname) const;

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FourBodyDecayConstructor & operator=(const FourBodyDecayConstructor &);

private:

  /**
   *  Whether or not to remove on-shell diagrams
   */
  unsigned int removeOnShell_;

  /**
   *  Option for the inclusion of intermediates
   */
  unsigned int interopt_;

  /**
   *  How to treat the widths of the intermediate particles
   */
  unsigned int widthopt_;

  /**
   * The minimum energy release for a three-body decay as a 
   * fraction of the parent mass
   */
  double minReleaseFraction_;

  /**
   *  Maximum number of EW gauge bosons
   */
  unsigned int maxBoson_;

  /**
   *  Maximum number of particles from the decaying particle list
   */
  unsigned int maxList_;

  /**
   *  Excluded Vertices
   */
  vector<VertexBasePtr> excludedVector_;

  /**
   *  Excluded Vertices
   */
  set<VertexBasePtr> excludedSet_;

};

}

#endif /* THEPEG_FourBodyDecayConstructor_H */
