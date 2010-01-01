// -*- C++ -*-
#ifndef HERWIG_VVHardGenerator_H
#define HERWIG_VVHardGenerator_H
//
// This is the declaration of the VVHardGenerator class.
//

#include "Herwig++/Shower/Base/HardestEmissionGenerator.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"
#include "Herwig++/MatrixElement/Powheg/VVKinematics.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "Herwig++/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {

using namespace ThePEG;
using namespace std;
/**
 * The VVHardGenerator class implements the hardest emission
 * in the POWHEG scheme for production of vector boson pairs.
 *
 * @see \ref VVHardGeneratorInterfaces "The interfaces"
 * defined for VVHardGenerator.
 */
class VVHardGenerator: public HardestEmissionGenerator {

  /**
   * Typedef for the BeamParticleData object
   */
  typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr;

public:

  /**
   * The default constructor.
   */
  VVHardGenerator();

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
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

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
   * @param pT The transverse momentum of the jet
   * @param yj The rapidity of the jet
   */
  double getResult(int emis_type, realVVKinematics R, Energy pT);
 
  /**
   *  generates the hardest emission (yj,p)
   * @param pnew The momenta of the new particles
   * @param emissiontype The type of emission, as for getResult
   * @return Whether not an emission was generated
   */
  bool getEvent(vector<Lorentz5Momentum> & pnew,unsigned int & emissiontype);
  
  /**
   *  sets the QCD, EW and PDF scales
   * @param pT The pT of the current step in the veto algorithm
   */
  void setTheScales(Energy pT);

  /**
   * The matrix element q + qb -> n + g times tk*uk 
   */
  Energy2 t_u_M_R_qqb_hel_amp(realVVKinematics R, bool getMatrix) const;


  /**
   * The matrix element q + g  -> n + q times tk*uk 
   */
  Energy2 t_u_M_R_qg_hel_amp(realVVKinematics R, bool getMatrix) const;

  /**
   * The matrix element g + qb -> n + q times tk*uk 
   */
  Energy2 t_u_M_R_gqb_hel_amp(realVVKinematics R, bool getMatrix) const;

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  double lo_me(bool getMatrix) const;

  /**
   * Recalculate hard vertex to include spin correlations for radiative events.
   */
  void recalculateVertex();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<VVHardGenerator> initVVHardGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VVHardGenerator & operator=(const VVHardGenerator &);

private:

  /**
   * If this boolean is true the n+1 body helicity amplitudes will be
   * used to calculate a hard vertex based on those kinematics for spin
   * correlations in the decays.
   */
  bool realMESpinCorrelations_;

  /**
   * Born / virtual 2->2 kinematics.
   */
  bornVVKinematics B_;

  /**
   * The colour & spin averaged n-body (leading order) matrix element squared.
   */
  double lo_me_;

  /**
   * The resolved 2->3 real emission kinematics.
   */
  realVVKinematics R_;

  /**
   * This specifies the emitting configuration: 
   * 1: q + qbar -> V1 + V2 + g
   * 2: q + g    -> V1 + V2 + q
   * 3: g + qbar -> V1 + V2 + qbar.
   */
  unsigned int channel_;

  /**
   * Identifies the space-like mother of the branching
   * as quark (+1) or antiquark (-1):
   */
  int fermionNumberOfMother_;

  /**
   * The radiative variable \tilde{x}.
   */
  double xt_;

  /**
   * The lower bound on the x integration \bar{x}. 
   */
  double xbar_;

  /**
   * The radiative variable y.
   */
  double y_;

  /**
   * The radiative variable theta_{2}.
   */
  double theta2_;

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr alphaS_;

  /**
   *  Constants for the sampling. The distribution is assumed to have the
   *  form \f$\frac{c}{{\rm GeV}}\times\left(\frac{{\rm GeV}}{p_T}\right)^n\f$ 
   */
  //@{
  /**
   * The power, \f$n\f$, for the sampling
   */
  double power_;

  /**
   *  The prefactor, \f$c\f$ for the \f$q\bar{q}\f$ channel
   */
  double preqqbar_;

  /**
   *  The prefactor, \f$c\f$ for the \f$qg\f$ channel
   */
  double preqg_;

  /**
   *  The prefactor, \f$c\f$ for the \f$g\bar{q}\f$ channel
   */
  double pregqbar_;

  /**
   * The QCD beta function divided by 4pi, (11-2/3*nf)/4/pi, with nf = 5.
   */
  double b0_;

  /**
   * The fundamental QCD scale in the one-loop alpha_{S} used for the crude
   * (not the very crude) overestimate of the Sudakov exponent. The default
   * value is set so such that alphaS(MZ), neglecting all flavour threshold
   * effects i.e. MZ*exp(-1/2/b0_/alphaS(MZ)).
   */
  Energy LambdaQCD_;

  /**
   *  The prefactors as a vector for easy use
   */
  vector<double> prefactor_;
  //@}

  /**
   *  Properties of the incoming particles
   */
  //@{
  /**
   *  Pointers to the ShowerProgenitor objects for the partons
   */
  ShowerProgenitorPtr qProgenitor_;
  ShowerProgenitorPtr qbProgenitor_;

  /**
   *  Pointers to the Shower particle objects for the partons
   */
  ShowerParticlePtr quark_;
  ShowerParticlePtr antiquark_;

  /**
   *  Pointers to the BeamParticleData objects
   */
  tcBeamPtr qHadron_;
  tcBeamPtr qbHadron_;
  //@}

  /**
   *  Properties of the boson and jets
   */
  //@{

  /**
   *  Pointers to the Shower particle objects for the partons
   */
  PPtr gluon_;
  PPtr V1_;
  PPtr V2_;
  PPtr emitted_;
  PPtr spacelikeSon_;
  vector<PPtr> children_;
  vector<PPtr> photons_;

  /**
   *  Flag indicating if the q & qbar are flipped or not i.e. this
   *  is true if q enters from the -z direction in the lab frame.
   */
  bool flipped_;

  /**
   *  the rapidity of the jet
   */
  double Yk_;

  /**
   *  The transverse momentum of the jet
   */
  Energy pT_;
  //@}

  /**
   *  The transverse momentum of the jet
   */
  Energy min_pT_;

  /**
   *  Option to impose helicity conservation on the real NLO ME's (to improve evaluation time).
   */
  bool helicityConservation_;

  // Work out the scales we want to use in the matrix elements and the pdfs:
  /**
   * Scale for alpha_S: pT^2 of the diboson system.
   */
   Energy2 QCDScale_;

  /**
   * Scale for real emission PDF: 
   */
  Energy2 PDFScale_;

  /**
   * Scale of electroweak vertices: mVV^2 the invariant mass of the diboson system.
   */
  Energy2 EWScale_;

  /**
   *  The vertices
   */
  AbstractFFVVertexPtr FFPvertex_;
  AbstractFFVVertexPtr FFWvertex_;
  AbstractFFVVertexPtr FFZvertex_;
  AbstractVVVVertexPtr WWWvertex_;
  AbstractFFVVertexPtr FFGvertex_;

  /**
   * A matrix element to hold information on the q + qbar -> V1 + V2 + g process
   */
  mutable ProductionMatrixElement qqb_hel_amps_;

  /**
   * A matrix element to hold information on the q + g    -> V1 + V2 + q process
   */
  mutable ProductionMatrixElement qg_hel_amps_;

  /**
   * A matrix element to hold information on the g + qbar -> V1 + V2 + qbar process
   */
  mutable ProductionMatrixElement gqb_hel_amps_;

  /**
   * A matrix element to hold information on the q + qbar -> V1 + V2 process
   */
  mutable ProductionMatrixElement lo_hel_amps_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of VVHardGenerator. */
template <>
struct BaseClassTrait<Herwig::VVHardGenerator,1> {
  /** Typedef of the first base class of VVHardGenerator. */
  typedef Herwig::HardestEmissionGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VVHardGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VVHardGenerator>
  : public ClassTraitsBase<Herwig::VVHardGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VVHardGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VVHardGenerator is implemented. It may also include several, space-separated,
   * libraries if the class VVHardGenerator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so HwPowhegMEHadron.so HwPowhegShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_VVHardGenerator_H */
