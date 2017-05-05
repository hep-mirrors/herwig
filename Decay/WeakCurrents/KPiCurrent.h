// -*- C++ -*-
//
// KPiCurrent.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_KPiCurrent_H
#define HERWIG_KPiCurrent_H
//
// This is the declaration of the KPiCurrent class.
//

#include "WeakDecayCurrent.h"
#include "Herwig/Utilities/Kinematics.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the KPiCurrent class.
 *
 * @see \ref KPiCurrentInterfaces "The interfaces"
 * defined for KPiCurrent.
 */
class KPiCurrent: public WeakDecayCurrent {

public:

  /**
   * The default constructor.
   */
  KPiCurrent();


  /** @name Methods for the construction of the phase space integrator. */
  //@{ 
  /**
   * Complete the construction of the decay mode for integration.
   * This version just adds the intermediate resonances, two outgoing mesons
   * and photon.
   * @param icharge The total charge of the outgoing particles in the current.
   * @param imode   The mode in the current being asked for.
   * @param mode    The phase space mode for the integration
   * @param iloc    The location of the of the first particle from the current in
   *                the list of outgoing particles.
   * @param ires    The location of the first intermediate for the current.
   * @param phase   The prototype phase space channel for the integration.
   * @param upp     The maximum possible mass the particles in the current are
   *                allowed to have.
   * @return Whether the current was sucessfully constructed.
   */
  virtual bool createMode(int icharge,unsigned int imode,DecayPhaseSpaceModePtr mode,
			  unsigned int iloc,unsigned int ires,
			  DecayPhaseSpaceChannelPtr phase,Energy upp);

  /**
   * The particles produced by the current. This just returns the two pseudoscalar
   * mesons and the photon.
   * @param icharge The total charge of the particles in the current.
   * @param imode The mode for which the particles are being requested
   * @param iq The PDG code for the quark
   * @param ia The PDG code for the antiquark
   * @return The external particles for the current.
   */
  virtual tPDVector particles(int icharge, unsigned int imode, int iq, int ia);
  //@}

  /**
   * Hadronic current. This version returns the hadronic current described above.
   * @param imode The mode
   * @param ichan The phase-space channel the current is needed for
   * @param scale The invariant mass of the particles in the current.
   * @param decay The decay products
   * @param meopt Option for the calculation of the matrix element
   * @return The current. 
   */
  virtual vector<LorentzPolarizationVectorE>  
  current(const int imode,const int ichan,Energy & scale,  
	  const ParticleVector & decay, DecayIntegrator::MEOption meopt) const;

  /**
   * Accept the decay. Checks the particles are the allowed mode.
   * @param id The id's of the particles in the current.
   * @return Can this current have the external particles specified.
   */
  virtual bool accept(vector<int> id);

  /**
   * Return the decay mode number for a given set of particles in the current. 
   * @param id The id's of the particles in the current.
   * @return The number of the mode
   */
  virtual unsigned int decayMode(vector<int> id);

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;

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

  /**
   * Breit-Wigner distributions
   */
  //@{
  /**
   * s-wave Breit-Wigner for the scalar resonances
   * @param q2 The scale
   * @param ires The resonances
   */
  Complex sWaveBreitWigner(Energy2 q2,unsigned int ires) const {
    Energy q=sqrt(q2),gam(ZERO);
    Energy2 m2=sqr(_scamass[ires]);
    if(q>_mK+_mpi) {
      Energy pX=Kinematics::pstarTwoBodyDecay(_scamass[ires],_mK,_mpi);
      Energy p =Kinematics::pstarTwoBodyDecay( q            ,_mK,_mpi);
      gam = _scawidth[ires]*m2/q2*p/pX;
    }
    return m2/(m2-q2-Complex(0.,1.)*q*gam);
  }
  
  /**
   *  p-wave Breit-Wigner for the vector resonances
   * @param q2 The scale
   * @param ires The resonances
   */
  Complex pWaveBreitWigner(Energy2 q2,unsigned int ires) const {
    Energy q=sqrt(q2),gam(ZERO);
    Energy2 m2=sqr(_vecmass[ires]);
    if(q>_mK+_mpi) {
      Energy pX=Kinematics::pstarTwoBodyDecay(_vecmass[ires],_mK,_mpi);
      Energy p =Kinematics::pstarTwoBodyDecay( q            ,_mK,_mpi);
      double ratio=p/pX;
      gam = _vecwidth[ires]*m2/q2*ratio*sqr(ratio);
    }
    return m2/(m2-q2-Complex(0.,1.)*q*gam);
  }
  //@}

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
   * Initialize this object after the setup phase before saving and
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
  static ClassDescription<KPiCurrent> initKPiCurrent;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  KPiCurrent & operator=(const KPiCurrent &);

private:

  /**
   *  Use local value of the parameters not those from the ParticleData objects
   */
  bool _localparameters;

  /**
   *  Whether to use \f$m^2\f$ or \f$Q^2\f$ in the projection operator.
   */
  bool _transverse;

  /**
   *  Normalizations of the vector and scalar pieces
   */
  //@{
  /**
   * \f$c_V\f$, normalization of the vector piece.
   */
  double _cV;

  /**
   * \f$c_S\f$, normalization of the scalar piece
   */
  double _cS;
  //@}

  /**
   * Parameters for the vector resonances
   */
  //@{
  /**
   *  Magnitude of the vector weights
   */
  vector<double> _vecmag;

  /**
   *  Phase of the vector weights
   */
  vector<double> _vecphase;

  /**
   *  Weights for the vector resonaces
   */
  vector<Complex> _vecwgt;

  /**
   *  Masses of the vector resonances
   */
  vector<Energy> _vecmass;

  /**
   *  Widths of the vector resonances
   */
  vector<Energy> _vecwidth;
  //@}

  /**
   * Parameters for the scalar resonances
   */
  //@{
  /**
   *  Magnitude of the scalar weights
   */
  vector<double> _scamag;

  /**
   *  Phase of the scalar weights
   */
  vector<double> _scaphase;

  /**
   *  Weights for the scalar resonances
   */
  vector<Complex> _scawgt;

  /**
   *  Masses of the scalar resonances
   */
  vector<Energy> _scamass;

  /**
   *  Widths of the scalar resonances
   */
  vector<Energy> _scawidth;
  //@}

  /**
   *  Masses for calculating the running widths
   */
  //@{
  /**
   * The pion mass
   */
  Energy _mpi;

  /**
   * The kaon mass
   */
  Energy _mK;
  //@}

  /**
   *  Map for the resonances
   */
  vector<int> _resmap;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of KPiCurrent. */
template <>
struct BaseClassTrait<Herwig::KPiCurrent,1> {
  /** Typedef of the first base class of KPiCurrent. */
  typedef Herwig::WeakDecayCurrent NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the KPiCurrent class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::KPiCurrent>
  : public ClassTraitsBase<Herwig::KPiCurrent> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::KPiCurrent"; }
  /**
   * The name of a file containing the dynamic library where the class
   * KPiCurrent is implemented. It may also include several, space-separated,
   * libraries if the class KPiCurrent depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwWeakCurrents.so"; }
};

/** @endcond */

}

#endif /* HERWIG_KPiCurrent_H */
