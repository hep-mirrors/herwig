// -*- C++ -*-
#ifndef Herwig_DynamicGluonMassGenerator_H
#define Herwig_DynamicGluonMassGenerator_H
//
// This is the declaration of the DynamicGluonMassGenerator class.
//



#include "Herwig/Hadronization/GluonMassGenerator.h"
#include "ClusterHadronizationHandler.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Cluster.h"


namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the GluonMassGenerator class.
 *
 * @see \ref GluonMassGeneratorInterfaces "The interfaces"
 * defined for GluonMassGenerator.
 */
class DynamicGluonMassGenerator: public GluonMassGenerator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DynamicGluonMassGenerator();

  /**
   * The destructor.
   */
  virtual ~DynamicGluonMassGenerator();
  //@}

public:


  /**
   * Returns the scale Qgtilde of the Sudakov in the gluon splitting
   */
  virtual Energy Qgtilde() const {
    return _Qgtilde; 
  }

  /**
   * Returns the minimum mass of a gluon that can by generated
   */
  virtual Energy minGluonMass() const { //ToDo: get it from hadron spectrum
    Energy m0, mu, md, ms;
    mu=getParticleData(ThePEG::ParticleID::u)->constituentMass();
    md=getParticleData(ThePEG::ParticleID::d)->constituentMass();
    ms=getParticleData(ThePEG::ParticleID::s)->constituentMass();
    m0=md;
    if(mu<m0){m0=mu;}
    if(ms<m0){m0=ms;}
    return m0;
  }


  /**
   * Generate a gluon mass drawn from the proposal distribution (overestimate)
   */
  virtual Energy generateProposal(Energy Qtilde) const;

  /**
   * The proposal distributions for the gluon mass (overestimate)
   */
  virtual InvEnergy PmgProposal(Energy mg, Energy mq) const;
  virtual InvEnergy PmgProposal(Energy mg) const;


  /**
   * The gluon mass distribution
   */
  virtual InvEnergy Pmg(Energy mg, Energy mq, Energy Qtilde) const;
  virtual InvEnergy Pmg(Energy mg, Energy Qtilde) const;

  
  /**
   * Generate a single gluon mass
   */
  virtual Energy generate(Energy Qtilde, Energy mgmax) const;
  virtual Energy generate(Energy Qtilde) const;
  virtual Energy generate() const;


  /**
   * Return the the strong coupling used in ClusterFission
   */
  double ClusterAlphaS(Energy2 q2) const { 
    if(q2>sqr(_clusteralphasfreeze)){return ClusterHadronizationHandler::currentHandler()->SM().alphaS(q2);}
    else{return ClusterHadronizationHandler::currentHandler()->SM().alphaS(sqr(_clusteralphasfreeze));}
   }

  

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:



  /**
   * Scale of the Sudakov in the gluon splitting
   */
  Energy _Qgtilde = 4.0*GeV;

  /**
   * freezing scale for the non-pert. alphas
   */
  Energy _clusteralphasfreeze = 1.0*GeV;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DynamicGluonMassGenerator & operator=(const DynamicGluonMassGenerator &);

};

}

#endif /* Herwig_DynamicGluonMassGenerator_H */
