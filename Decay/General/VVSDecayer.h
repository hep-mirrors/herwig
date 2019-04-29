// -*- C++ -*-
#ifndef THEPEG_VVSDecayer_H
#define THEPEG_VVSDecayer_H
//
// This is the declaration of the VVSDecayer class.
//

#include "GeneralTwoBodyDecayer.h"
#include "ThePEG/Helicity/Vertex/Scalar/VVSVertex.h"
#include "ThePEG/Repository/EventGenerator.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::VVSVertexPtr;

/** \ingroup Decay
 * The VVSDecayer class implements the decay of a vector to a
 * vector and a scalar in a general model. It holds an VVSVertex pointer
 * that must be typecast from the  VertexBase pointer helid in the
 * GeneralTwoBodyDecayer. It implents the virtual functions me2() and
 * partialWidth(). 
 *
 * @see \ref VVSDecayerInterfaces "The interfaces"
 * defined for VVSDecayer.
 */
class VVSDecayer: public GeneralTwoBodyDecayer {

public:

  /**
   * The default constructor.
   */
  VVSDecayer() {}

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Return the matrix element squared for a given mode and phase-space channel
   * @param ichan The channel we are calculating the matrix element for.
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(const int ichan, const Particle & part,
		     const ParticleVector & decay, MEOption meopt) const;

  /**
   * Function to return partial Width
   * @param inpart The decaying particle.
   * @param outa One of the decay products.
   * @param outb The other decay product.
   */
  virtual Energy partialWidth(PMPair inpart, PMPair outa, 
			      PMPair outb) const;

  /**
   *  Set the information on the decay
   */
  virtual void setDecayInfo(PDPtr incoming, PDPair outgoing,
			    vector<VertexBasePtr>,
			    map<ShowerInteraction,VertexBasePtr> &,
			    const vector<map<ShowerInteraction,VertexBasePtr> > &,
			    map<ShowerInteraction,VertexBasePtr>);
  
  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection()  {
    POWHEGType output = FSR;
    for(auto vertex : vertex_) {
      if(vertex->orderInAllCouplings()!=1) {
	output = No;
	break;
      }
    }
    return output;
  }

  /**
   *  Three-body matrix element including additional QCD radiation
   */
  virtual double threeBodyME(const int , const Particle & inpart,
			     const ParticleVector & decay,
			     ShowerInteraction inter, MEOption meopt);
  //@}

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
  VVSDecayer & operator=(const VVSDecayer &) = delete;

private:

  /**
   *  Abstract pointer to AbstractVVSVertex
   */
  vector<AbstractVVSVertexPtr> vertex_;

  /**
   * Pointer to the perturbative vertex
   */
  vector<VVSVertexPtr> perturbativeVertex_;

  /**
   *  Abstract pointer to AbstractVVVVertex for QCD radiation from incoming vector
   */
  map<ShowerInteraction,AbstractVVVVertexPtr> incomingVertex_;

  /**
   *  Abstract pointer to AbstractVVVVertex for QCD radiation from the first outgoing vector
   */
  map<ShowerInteraction,AbstractVVVVertexPtr> outgoingVertexV_;

  /**
   *  Abstract pointer to AbstractVVVVertex for QCD radiation from the second outgoing vector
   */
  map<ShowerInteraction,AbstractVSSVertexPtr> outgoingVertexS_;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix rho_;

  /**
   *  Vector wavefunctions
   */
  mutable vector<Helicity::VectorWaveFunction> vectors_[2];
  
private:

  /**
   *  Members for the POWHEG correction
   */
  //@{
  /**
   *  Spin density matrix for 3 body decay
   */
  mutable RhoDMatrix rho3_;

  /**
   *  Vector wavefunctions
   */
  mutable vector<Helicity::VectorWaveFunction> vectors3_[2];

    /**
   *  Vector wavefunction for 3 body decay
   */
  mutable vector<Helicity::VectorWaveFunction> gluon_;
  //@}
};

}

#endif /* THEPEG_VVSDecayer_H */
