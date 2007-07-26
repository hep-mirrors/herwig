// -*- C++ -*-
#ifndef HERWIG_DtoPiPiPiFOCUS_H
#define HERWIG_DtoPiPiPiFOCUS_H
//
// This is the declaration of the DtoPiPiPiFOCUS class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "DtoPiPiPiFOCUS.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the DtoPiPiPiFOCUS class.
 *
 * @see \ref DtoPiPiPiFOCUSInterfaces "The interfaces"
 * defined for DtoPiPiPiFOCUS.
 */
class DtoPiPiPiFOCUS: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  DtoPiPiPiFOCUS();

  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param dm The decay mode
   */
  virtual int modeNumber(bool & cc,const DecayMode & dm) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param vertex Output the information on the vertex for spin correlations
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2(bool vertex, const int ichan,const Particle & part,
	     const ParticleVector & decay) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

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


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DtoPiPiPiFOCUS> initDtoPiPiPiFOCUS;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DtoPiPiPiFOCUS & operator=(const DtoPiPiPiFOCUS &);

private:

  /**
   *  Parameters for the K-matrix
   */
  //@{
  /**
   *  Masses of the poles for the K-matrix
   */
  vector<Energy> _malpha;

  /**
   *  The \f$g_{\pi\pi}\f$ coupling
   */
  vector<Energy> _gpipi;

  /**
   *  The \f$g_{K\bar{K}}\f$ coupling
   */
  vector<Energy> _gKK;

  /**
   *  The \f$g_{4\pi}\f$ coupling
   */
  vector<Energy> _g4pi;

  /*
   *  The \f$g_{\eta\eta}\f$ coupling
   */
  vector<Energy> _getaeta;

  /**
   *  The \f$g_{\eta\eta'}\f$ coupling
   */
  vector<Energy> _getaetap;

  /**
   *  The g couplings for easy access
   */
  vector<vector<Energy> > _gcoup;

  /**
   *  \f$s^{\rm scatt}_0\f$
   */
  Energy2 _s0scatt;

  /**
   *  \f$s_{A_0}\f$
   */
  Energy2 _sA0;

  /**
   * \f$s_A\f$
   */
  double _sA;

  /**
   *  Pion mass
   */
  Energy _mpi;

  /**
   *  Kaon mass
   */
  Energy _mK;

  /**
   * \f$\eta\f$ mass
   */
  Energy _meta;

  /**
   *  \f$\eta'\f$ mass
   */
  Energy _metap;

  /**
   *  The magnitude of the \f$f_{\rm prod}\f$ couplings for \f$D\f$
   */
  vector<double> _fprodmagD;

  /**
   *  The phase of the \f$f_{\rm prod}\f$ couplings for \f$D\f$
   */
  vector<double> _fprodphaseD;

  /**
   * The \f$f_{\rm prod}\f$ couplings for \f$D\f$
   */
  vector<Complex> _fprodD;

  /**
   *  The magnitude of the \f$f_{\rm prod}\f$ couplings for \f$D_s\f$
   */
  vector<double> _fprodmagDs;

  /**
   *  The phase of the \f$f_{\rm prod}\f$ couplings for \f$D_s\f$
   */
  vector<double> _fprodphaseDs;

  /**
   * The \f$f_{\rm prod}\f$ couplings for \f$D_s\f$
   */
  vector<Complex> _fprodDs;

  /**
   *  The magnitude of the \f$\beta\f$ couplings for \f$D\f$
   */
  vector<Energy> _betamagD;

  /**
   *  The phase of the \f$\beta\f$ couplings for \f$D\f$
   */
  vector<double> _betaphaseD;

  /**
   * The \f$\beta\f$ couplings for \f$D\f$
   */
  vector<complex<Energy> > _betaD;

  /**
   *  The magnitude of the \f$\beta\f$ couplings for \f$D_s\f$
   */
  vector<Energy> _betamagDs;

  /**
   *  The phase of the \f$\beta\f$ couplings for \f$D_s\f$
   */
  vector<double> _betaphaseDs;

  /**
   * The \f$\beta\f$ couplings for \f$D_s\f$
   */
  vector<complex<Energy> > _betaDs;
  //@}

  /**
   *  Amplitudes and phases for the vector and tensor mesons
   */
  //@{
  /**
   *  Magnitude for the \f$f_2(1270)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
   */
  InvEnergy2 _aDsf2;

  /**
   *  Phase for the \f$f_2(1270)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
   */
  double _phiDsf2;

  /**
   *  Magnitude for the \f$\rho^0(1450)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
   */
  double _aDsrho1450;

  /**
   *  Phase for the \f$\rho^0(1450)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
   */
  double _phiDsrho1450;

  /**
   *  Amplitude for the \f$f_2(1270)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
   */
  complex<InvEnergy2> _cDsf2;

  /**
   *  Amplitude for the \f$\rho^0(1450)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
   */
  Complex _cDsrho1450;

  /**
   *  Magnitude for the \f$f_2(1270)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
   */
  InvEnergy2 _aDf2;

  /**
   *  Phase for the \f$f_2(1270)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
   */
  double _phiDf2;

  /**
   *  Magnitude for the \f$\rho^0(770)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
   */
  double _aDrho770;

  /**
   *  Phase for the \f$\rho^0(770)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
   */
  double _phiDrho770;

  /**
   *  Amplitude for the \f$f_2(1270)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
   */
  complex<InvEnergy2> _cDf2;

  /**
   *  Amplitude for the \f$\rho^0(770)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
   */
  Complex _cDrho770;
  //@}


  /**
   *  Masses and Widths for the vector and tensor mesons
   */
  //@{
  /**
   *  Mass of the \f$\rho(770)\f$
   */
  Energy _mrho770;

  /**
   *  Mass of the \f$\rho(1450)\f$
   */
  Energy _mrho1450;

  /**
   *  Mass of the \f$f_2(1270)\f$
   */
  Energy _mf2;

  /**
   *  Width of the \f$\rho(770)\f$
   */
  Energy _wrho770;

  /**
   *  Width of the \f$\rho(1450)\f$
   */
  Energy _wrho1450;

  /**
   *  Width of the \f$f_2(1270)\f$
   */
  Energy _wf2;
  //@}

  /**
   * Parameters for the phase-space integration
   */
  //@{
  /**
   *  Maximum weights for the different modes
   */
  vector<double> _maxwgt;

  /**
   *  Weights for the different phase-space channels
   */
  vector<double> _weights;
  //@}

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DtoPiPiPiFOCUS. */
template <>
struct BaseClassTrait<Herwig::DtoPiPiPiFOCUS,1> {
  /** Typedef of the first base class of DtoPiPiPiFOCUS. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DtoPiPiPiFOCUS class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DtoPiPiPiFOCUS>
  : public ClassTraitsBase<Herwig::DtoPiPiPiFOCUS> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DtoPiPiPiFOCUS"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DtoPiPiPiFOCUS is implemented. It may also include several, space-separated,
   * libraries if the class DtoPiPiPiFOCUS depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSMDecay.so"; }
};

/** @endcond */

}

#include "DtoPiPiPiFOCUS.icc"

#endif /* HERWIG_DtoPiPiPiFOCUS_H */
