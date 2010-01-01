// -*- C++ -*-
#ifndef HERWIG_DrellYanHardGenerator_H
#define HERWIG_DrellYanHardGenerator_H
//
// This is the declaration of the DrellYanHardGenerator class.
//

#include "Herwig++/Shower/Base/HardestEmissionGenerator.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"

namespace Herwig {

using namespace ThePEG;
using namespace std;
/**
 * The DrellYanHardGenerator class implements the hardest emission
 * in the POWHEG scheme for Drell-Yan production of vector bosons.
 *
 * @see \ref DrellYanHardGeneratorInterfaces "The interfaces"
 * defined for DrellYanHardGenerator.
 */
class DrellYanHardGenerator: public HardestEmissionGenerator {

  /**
   * Typedef for the BeamParticleData object
   */
  typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr;

public:

  /**
   * The default constructor.
   */
  DrellYanHardGenerator();

  /**
   *  Implementation of virtual members from HardestEmissionGenerator
   */
  //@{
  /**
   *  Member to generate the hardest emission
   */
  virtual HardTreePtr generateHardest(ShowerTreePtr);

  /**
   *  Member to decide if the inheriting class can handle this process
   */
  virtual bool canHandle(ShowerTreePtr);
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
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

protected:

  /**
   * Returns the matrix element for a given type of process,
   * rapidity of the jet \f$y_j\f$ and transverse momentum \f$p_T\f$
   * @param emis_type the type of emission,
   * (0 is \f$q\bar{q}\to Vg\f$, 1 is \f$qg\to Vq\f$ and 2 is \f$g\bar{q}\to V\bar{q}\f$)
   * @param pt The transverse momentum of the jet
   * @param yj The rapidity of the jet
   */
  double getResult(int emis_type, Energy pt, double yj);
 
  /**
   *  generates the hardest emission (yj,p)
   * @param pnew The momenta of the new particles
   * @param emissiontype The type of emission, as for getResult
   * @return Whether not an emission was generated
   */
  bool getEvent(vector<Lorentz5Momentum> & pnew,int & emissiontype);
  
private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DrellYanHardGenerator> initDrellYanHardGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DrellYanHardGenerator & operator=(const DrellYanHardGenerator &);

private:

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr _alphaS;

  /**
   *  Constants for the sampling. The distribution is assumed to have the
   *  form \f$\frac{c}{{\rm GeV}}\times\left(\frac{{\rm GeV}}{p_T}\right)^n\f$ 
   */
  //@{
  /**
   * The power, \f$n\f$, for the sampling
   */
  double _power;

  /**
   *  The prefactor, \f$c\f$ for the \f$q\bar{q}\f$ channel
   */
  double _preqqbar;

  /**
   *  The prefactor, \f$c\f$ for the \f$qg\f$ channel
   */
  double _preqg;

  /**
   *  The prefactor, \f$c\f$ for the \f$g\bar{q}\f$ channel
   */
  double _pregqbar;

  /**
   *  The prefactors as a vector for easy use
   */
  vector<double> _prefactor;
  //@}

  /**
   *  Properties of the incoming particles
   */
  //@{
  /**
   *  Pointers to the BeamParticleData objects
   */
  vector<tcBeamPtr> _beams;
  
  /**
   *  Pointers to the ParticleDataObjects for the partons
   */
  vector<tcPDPtr> _partons;
  //@}

  /**
   *  Properties of the boson and jets
   */
  //@{
  /**
   *  The rapidity of the gauge boson
   */
  double _yb;

  /**
   *  The mass of the gauge boson
   */
  Energy _mass;

  /**
   *  Whether the quark is in the + or - z direction
   */
  bool _quarkplus;

  /**
   *  the rapidity of the jet
   */
  double _yj;

  /**
   *  The transverse momentum of the jet
   */
  Energy _pt;
  //@}

  /**
   *  The transverse momentum of the jet
   */
  Energy _min_pt;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DrellYanHardGenerator. */
template <>
struct BaseClassTrait<Herwig::DrellYanHardGenerator,1> {
  /** Typedef of the first base class of DrellYanHardGenerator. */
  typedef Herwig::HardestEmissionGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DrellYanHardGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DrellYanHardGenerator>
  : public ClassTraitsBase<Herwig::DrellYanHardGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DrellYanHardGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DrellYanHardGenerator is implemented. It may also include several, space-separated,
   * libraries if the class DrellYanHardGenerator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so HwPowhegMEHadron.so HwPowhegShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DrellYanHardGenerator_H */
