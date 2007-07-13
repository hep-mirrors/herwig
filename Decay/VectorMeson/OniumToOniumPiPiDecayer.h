// -*- C++ -*-
#ifndef HERWIG_OniumToOniumPiPiDecayer_H
#define HERWIG_OniumToOniumPiPiDecayer_H
//
// This is the declaration of the OniumToOniumPiPiDecayer class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "OniumToOniumPiPiDecayer.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Decay
 *  The OniumToOniumPiPiDecayer class uses the matrix element of 
 *  Brown and Cahn PRL35, 1 (1975) for the decay of onium resonaces to
 *  lighter states and pion pairs. The matrix element is given by
 * \f[\mathcal{M} = \epsilon'\cdot\epsilon\left[
 *    \mathcal{A}\left(q^2-2m^2_\pi\right)+\mathcal{B}E_1E_2\right]
 *    +\mathcal{C}\left((\epsilon'\cdot q_1)(\epsilon\cdot q_2)+
 *                     (\epsilon'\cdot q_2)(\epsilon\cdot q_1)\right),\f]
 * where \f$\epsilon'\f$ is the polarization vector of the decaying onium resonance,
 *       \f$\epsilon\f$  is the polarization vector of the outgoing onium resonance,
 *    \f$\mathcal{A}\f$, \f$\mathcal{B}\f$ and \f$\mathcal{C}\f$ are complex couplings,
 *    \f$m_\pi\f$ is the pion mass, \f$E_{1,2}\f$ are the pion energies, \f$q_{1,2}\f$ 
 *    are the pion momenta and \f$q\f$ is the momentum of the \f$\pi\pi\f$ system.
 * 
 * The results of hep-ex/9909038 are used for \f$\psi'\to J/\psi\f$ and arXiv:0706.2317
 * for \f$\Upsilon(3S)\f$ and \f$\Upsilon(2S)\f$ decays.
 * The remaining parameters are choosen
 * to approximately reproduce the distributions from hep-ex/0604031 and hep-ex/0508023.
 *
 * @see \ref OniumToOniumPiPiDecayerInterfaces "The interfaces"
 * defined for OniumToOniumPiPiDecayer.
 */
class OniumToOniumPiPiDecayer: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  OniumToOniumPiPiDecayer();
  
public:
  
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

  /**
   * Initialize this object to the begining of the run phase.
   */
  inline virtual void doinitrun();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<OniumToOniumPiPiDecayer> initOniumToOniumPiPiDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  OniumToOniumPiPiDecayer & operator=(const OniumToOniumPiPiDecayer &);

private:

  /**
   * the PDG codes for the incoming onium resonace
   */
  vector<long> _incoming;

  /**
   * the PDG codes for the outgoing onium resonance
   */
  vector<long> _outgoing;

  /**
   * the maximum weight for the integration
   */
  vector<double> _maxweight;

  /**
   *  Initial size of the vectors
   */
  unsigned int _initsize;

  /**
   *  The couplings for the decays
   */
  //@{
  /**
   *  The real part of \f$A\f$
   */
  vector<InvEnergy2> _reA;

  /**
   *  The imaginary part of \f$A\f$
   */
  vector<InvEnergy2> _imA;

  /**
   *  The complex \f$A\f$ coupling
   */
  vector<complex<InvEnergy2> > _cA;

  /**
   *  The real part of \f$B\f$
   */
  vector<InvEnergy2> _reB;

  /**
   *  The imaginary part of \f$B\f$
   */
  vector<InvEnergy2> _imB;

  /**
   *  The complex \f$B\f$ coupling
   */
  vector<complex<InvEnergy2> > _cB;

  /**
   *  The real part of \f$C\f$
   */
  vector<InvEnergy2> _reC;

  /**
   *  The imaginary part of \f$C\f$
   */
  vector<InvEnergy2> _imC;

  /**
   *  The complex \f$C\f$ coupling
   */
  vector<complex<InvEnergy2> > _cC;
  //@}

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of OniumToOniumPiPiDecayer. */
template <>
struct BaseClassTrait<Herwig::OniumToOniumPiPiDecayer,1> {
  /** Typedef of the first base class of OniumToOniumPiPiDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the OniumToOniumPiPiDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::OniumToOniumPiPiDecayer>
  : public ClassTraitsBase<Herwig::OniumToOniumPiPiDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::OniumToOniumPiPiDecayer"; }
  /**
   * The name of a file containing the dynamic library where the class
   * OniumToOniumPiPiDecayer is implemented. It may also include several, space-separated,
   * libraries if the class OniumToOniumPiPiDecayer depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwVMDecay.so"; }
};

/** @endcond */

}

#include "OniumToOniumPiPiDecayer.icc"

#endif /* HERWIG_OniumToOniumPiPiDecayer_H */
