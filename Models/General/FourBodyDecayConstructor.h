// -*- C++ -*-
#ifndef THEPEG_FourBodyDecayConstructor_H
#define THEPEG_FourBodyDecayConstructor_H
//
// This is the declaration of the FourBodyDecayConstructor class.
//

#include "NBodyDecayConstructorBase.h"
#include "ThePEG/Helicity/Vertex/VertexBase.h"
#include "Herwig/Decay/General/GeneralFourBodyDecayer.fh"
#include "PrototypeVertex.h"

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
    interOpt_(0), widthOpt_(1), particleType_(false) {}

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
  void createDecayMode(vector<NBDiagram> &,bool,double);

  /**
   * Create the decayer
   * @param diagrams The diagrams for the decay
   * @param inter Option for intermediates
   */
  GeneralFourBodyDecayerPtr createDecayer(vector<NBDiagram> & diagrams,
					  bool inter, double symfac) const;

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FourBodyDecayConstructor & operator=(const FourBodyDecayConstructor &);

private:

  /**
   *  Option for the inclusion of intermediates
   */
  unsigned int interOpt_;

  /**
   *  How to treat the widths of the intermediate particles
   */
  unsigned int widthOpt_;

  /**
   *  Particles to override the default list
   */
  vector<PDPtr> particles_;

  /**
   *  Types of particles
   */
  bool particleType_;
};

}

#endif /* THEPEG_FourBodyDecayConstructor_H */
