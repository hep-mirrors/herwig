// -*- C++ -*-
//
// DipoleSplittingKinematics.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleSplittingKinematics_H
#define HERWIG_DipoleSplittingKinematics_H
//
// This is the declaration of the DipoleSplittingKinematics class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"

#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Utilities/UtilityBase.h"

#include "Herwig/Shower/Dipole/Utility/DipoleMCCheck.h"

namespace Herwig {

  using namespace ThePEG;

  class DipoleIndex;
  class DipoleSplittingInfo;
  class DipoleSplittingKernel;

  /**
   * \ingroup DipoleShower
   * \author Simon Platzer
   *
   * \brief DipoleSplittingKinematics is the base class for dipole splittings
   * as performed in the dipole shower.
   *
   * @see \ref DipoleSplittingKinematicsInterfaces "The interfaces"
   * defined for DipoleSplittingKinematics.
   */
  class DipoleSplittingKinematics: public HandlerBase {

  public:

    /** @name Standard constructors and destructors. */
    //@{
    /**
     * The default constructor.
     */
    DipoleSplittingKinematics();

    /**
     * The destructor.
     */
    virtual ~DipoleSplittingKinematics();
    //@}

  public:

    /**
     * Return the boundaries in between the evolution
     * variable random number is to be sampled; the lower
     * cuoff is assumed to correspond to the infrared cutoff.
     */
    virtual pair<double,double> kappaSupport(const DipoleSplittingInfo&) const {
      return {0.0,1.0};
    }

    /**
     * Return the boundaries in between the momentum
     * fraction random number is to be sampled.
     */
    virtual pair<double,double> xiSupport(const DipoleSplittingInfo&) const {
      return {0.0,1.0};
    }

    /**
     * Return the dipole scale associated to the
     * given pair of emitter and spectator. This
     * should be the invariant mass or absolute value
     * final/final or initial/initial and the absolute
     * value of the momentum transfer for intial/final or
     * final/initial dipoles.
     */
    virtual Energy dipoleScale(const Lorentz5Momentum& pEmitter,
			       const Lorentz5Momentum& pSpectator) const {
	// MEMinBias produces non-zero zeros.
      if(abs(pEmitter*pSpectator)<0.0000001*GeV2)return ZERO;
      assert(pEmitter*pSpectator >= ZERO);
      return sqrt(2.*pEmitter*pSpectator);
    }

    /**
     * Return the mass of the system absorbing 
     * the recoil in the dipole splitting. 
     * This is overloaded in the decay dipoles.
     */
    virtual Energy recoilMassKin(const Lorentz5Momentum&,
			      const Lorentz5Momentum& pSpectator) const {
      return pSpectator.m();
    }

    /**
     * Return the maximum pt for the given dipole scale.
     */
    virtual Energy ptMax(Energy dScale, 
			 double emX, double specX,
			 const DipoleIndex& dIndex,
			 const DipoleSplittingKernel& split) const =0;

    /**
     * Return the maximum pt for the given dipole scale.
     */
    virtual Energy ptMax(Energy dScale, 
			 double emX, double specX,
			 const DipoleSplittingInfo& dInfo,
			 const DipoleSplittingKernel& split) const; 

    /**
     * Return the maximum virtuality for the given dipole scale.
     */
    virtual Energy QMax(Energy dScale, 
			double emX, double specX,
			const DipoleIndex& dIndex,
			const DipoleSplittingKernel& split) const =0;

    /**
     * Return the maximum virtuality for the given dipole scale.
     */
    virtual Energy QMax(Energy dScale, 
			double emX, double specX,
			const DipoleSplittingInfo& dInfo,
			const DipoleSplittingKernel& split) const;

    /**
     * Return the pt given a virtuality.
     */
    virtual Energy PtFromQ(Energy scale, const DipoleSplittingInfo&) const = 0;

    /**
     * Return the virtuality given a pt.
     */
    virtual Energy QFromPt(Energy scale, const DipoleSplittingInfo&) const = 0;

    /**
     * Return the infrared cutoff.
     */
    virtual Energy IRCutoff() const { return theIRCutoff; }

    /**
     * Return the minimum momentum fraction for
     * incoming partons
     */
    double xMin() const { return theXMin; }

    /**
     * Generate a pt
     */
    Energy generatePt(double r, Energy dScale,
		      double emX, double specX,
		      const DipoleIndex& dIndex,
		      const DipoleSplittingKernel& split,
		      double& weight) const;

    /**
     * Return the random number associated to
     * the given pt.
     */
    virtual double ptToRandom(Energy pt, Energy dScale,
			      double emX, double specX,
			      const DipoleIndex& dIndex,
			      const DipoleSplittingKernel& split) const;

    /**
     * Return the boundaries on the momentum fraction
     */
    virtual pair<double,double> zBoundaries(Energy pt,
					    const DipoleSplittingInfo& dInfo,
					    const DipoleSplittingKernel& split) const = 0;

    /**
     * Enumerate the variants of sampling z
     */
    enum ZSamplingOptions {

      FlatZ = 0,
      OneOverZ,
      OneOverOneMinusZ,
      OneOverZOneMinusZ

    };

    /**
     * Generate a z value flat
     */
    double generateZ(double r, Energy pt, int sampling,
		     const DipoleSplittingInfo& dInfo,
		     const DipoleSplittingKernel& split,
		     double& weight) const;

    /**
     * Generate splitting variables given three random numbers
     * and the momentum fractions of the emitter and spectator.
     * Return true on success.
     */
    virtual bool generateSplitting(double kappa, double xi, double phi,
				   DipoleSplittingInfo& info,
				   const DipoleSplittingKernel&) = 0;

    /**
     * Get the splitting phasespace weight associated to
     * the last call to generateSplitting. This is taken to
     * be the single particle phasespace times 16 \pi^2 divided
     * by the relevant propagator invariant.
     */
    double jacobian() const { return theJacobian; }

    /**
     * Return true, if this splitting kinematics
     * class is capable of delivering an overestimate
     * to the jacobian.
     */
    virtual bool haveOverestimate() const { return false; }

    /**
     * Return the overestimated jacobian for the
     * last generated parameters.
     */
    virtual double jacobianOverestimate() const { return -1.; }

    /**
     * Return the last generated pt
     */
    Energy lastPt() const { return theLastPt; }

    /**
     * Return the last generated momentum fraction.
     */
    double lastZ() const { return theLastZ; }

    /**
     * Return the last calculated zPrime for massive FF and decay dipoles.
     */
    
    // Do not need in current implementation,
    // using lastSplittingParameters instead.
    //double lastZPrime() const { return theLastZPrime; }

    /**
     * Return the last generated azimuthal angle.
     */
    double lastPhi() const { return theLastPhi; }

    /**
     * Return the momentum fraction, by which the emitter's
     * momentum fraction should be divided after the splitting.
     */
    double lastEmitterZ() const { return theLastEmitterZ; }

    /**
     * Return the momentum fraction, by which the spectator's
     * momentum fraction should be divided after the splitting.
     */
    double lastSpectatorZ() const { return theLastSpectatorZ; }

    /**
     * Return any additional parameters needed to
     * evaluate the splitting kernel or to generate the 
     * full splitting.
     */
    const vector<double>& lastSplittingParameters() const { return theLastSplittingParameters; }

    /**
     * Complete a DipoleSplittingInfo object with
     * the parameters generated by the last call to
     * generateSplitting()
     */
    void prepareSplitting(DipoleSplittingInfo& dInfo);

  public:

    /**
     * Generate the full kinematics given emitter and
     * spectator momentum and a previously completeted
     * DipoleSplittingInfo object.
     */
    virtual void generateKinematics(const Lorentz5Momentum& pEmitter,
				    const Lorentz5Momentum& pSpectator,
				    const DipoleSplittingInfo& dInfo) = 0;



    /**
     * Return the emitter's momentum after the splitting.
     */
    const Lorentz5Momentum& lastEmitterMomentum() const { return theEmitterMomentum; }

    /**
     * Return the spectator's momentum after the splitting.
     */
    const Lorentz5Momentum& lastSpectatorMomentum() const { return theSpectatorMomentum; }

    /**
     * Return the emission's momentum.
     */
    const Lorentz5Momentum& lastEmissionMomentum() const { return theEmissionMomentum; }

    /*
     * Return true, if there is a transformation which should
     * be applied to all other final state particles except the ones
     * involved in the splitting after having performed the splitting.
     */
    virtual bool doesTransform () const { return false; }

    /*
     * Use the Dipole scale instead of hardpt for z-boundaries.
     */
    int openZBoundaries() const { return theOpenZBoundaries; }

    /*
     * perform the transformation if required.
     */
    virtual Lorentz5Momentum transform (const Lorentz5Momentum& p) const { return p; }

    /*
     * Return true if this splitting is of a dipole which contains
     * a decayed parton and requires the remnant to absorb the recoil.
     */
    virtual bool isDecay() const { return false; }

    /**
     * Perform the recoil in the case of a decayed parton
     */
    //virtual Lorentz5Momentum decayRecoil ( const Lorentz5Momentum& p, const int) { return p; }

    /**
     * Perform the recoil in the case of a decayed parton
     */
    virtual void decayRecoil ( PList& ) {};

 // {;}

  protected:

    /**
     * Calculate a transverse momentum for the given momenta,
     * invariant pt and azimuth.
     */
    Lorentz5Momentum getKt(const Lorentz5Momentum& p1,
			   const Lorentz5Momentum& p2,
			   Energy pt,
			   double phi,
			   bool spacelike = false) const;

    /**
     * Set the splitting phasespace weight associated to
     * the last call to generateSplitting. This is taken to
     * be the single particle phasespace times 16 \pi^2 divided
     * by the relevant propagator invariant.
     */
    void jacobian(double w) { theJacobian = w; }

    /**
     * Set the last generated pt
     */
    void lastPt(Energy p) { theLastPt = p; }

    /**
     * Set the last generated momentum fraction.
     */
    void lastZ(double z) { theLastZ = z; }

    /**
     * Set the last calculated zPrime for massive FF and decay dipoles.
     */
    // Do not need in current implementation,
    // using lastSplittingParameters instead.
    //void lastZPrime(double zPrime) { theLastZPrime = zPrime; }

    /**
     * Set the last generated azimuthal angle.
     */
    void lastPhi(double p) { theLastPhi = p; }

    /**
     * Set the momentum fraction, by which the emitter's
     * momentum fraction should be divided after the splitting.
     */
    void lastEmitterZ(double z) { theLastEmitterZ = z; }

    /**
     * Set the momentum fraction, by which the spectator's
     * momentum fraction should be divided after the splitting.
     */
    void lastSpectatorZ(double z) { theLastSpectatorZ = z; }

    /**
     * Access any additional parameters needed to
     * evaluate the splitting kernel or to generate the 
     * full splitting.
     */
    vector<double>& splittingParameters() { return theLastSplittingParameters; }

    /**
     * Set the emitter's momentum after the splitting.
     */
    void emitterMomentum(const Lorentz5Momentum& p) { theEmitterMomentum = p; }

    /**
     * Set the spectator's momentum after the splitting.
     */
    void spectatorMomentum(const Lorentz5Momentum& p) { theSpectatorMomentum = p; }

    /**
     * Set the emission's momentum.
     */
    void emissionMomentum(const Lorentz5Momentum& p) { theEmissionMomentum = p; }

    /**
     * Set the momentum of the recoil system before the splitting.
     * 25/05/2016 - Not currently used.
     */
    //void recoilMomentum( const Lorentz5Momentum& mom ) { theRecoilMomentum = mom; }

    /**
     * Set the momentum of the recoil system after the splitting.
     */
    void splitRecoilMomentum( const Lorentz5Momentum& mom ) { theSplitRecoilMomentum = mom; }

    /**
     * Return the momentum of the recoil system before splitting.
     * 25/05/2016 - Not currently used.
     */
    //const Lorentz5Momentum& recoilMomentum() const { return theRecoilMomentum; }

    /**
     * Return the momentum of the recoil system after splitting.
     */
    const Lorentz5Momentum& splitRecoilMomentum() const { return theSplitRecoilMomentum; }

  public:

    /**
     * The standard Init function used to initialize the interfaces.
     * Called exactly once for each class by the class description system
     * before the main function starts or
     * when this class is dynamically loaded.
     */
    static void Init();

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


    // If needed, insert declarations of virtual function defined in the
    // InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

  private:

    /**
     * The infrared cutoff associated to this
     * splitting kinematics.
     */
    Energy theIRCutoff;

    /**
     * The minimum momentum fraction for
     * incoming partons
     */
    double theXMin;

    /**
     * The last calculated splitting phase space weight.
     */
    double theJacobian;

    /**
     * The last generated pt
     */
    Energy theLastPt;

    /**
     * The last generated momentum fraction.
     */
    double theLastZ;

    /**
     * The last calculated zPrime required for massive FF
     * and decay kinematics dipoles.
     * zPrime := qi.nk / (qi+qj).nk (qj = emission momentum)
     */
    // Do not need in current implementation,
    // using lastSplittingParameters instead.
    //double theLastZPrime;

    /**
     * The last generated azimuthal angle.
     */
    double theLastPhi;

    /**
     * The momentum fraction, by which the emitter's
     * momentum fraction should be divided after the splitting.
     */
    double theLastEmitterZ;

    /**
     * The momentum fraction, by which the spectator's
     * momentum fraction should be divided after the splitting.
     */
    double theLastSpectatorZ;

    /**
     * Any additional parameters needed to
     * evaluate the splitting kernel or to generate the 
     * full splitting.
     */
    vector<double> theLastSplittingParameters;

    /**
     * The emitter's momentum after the splitting.
     */
    Lorentz5Momentum theEmitterMomentum;

    /**
     * The emission's momentum after the splitting.
     */
    Lorentz5Momentum theEmissionMomentum;

    /**
     * The spectator's momentum after the splitting.
     */
    Lorentz5Momentum theSpectatorMomentum;

    /**
     * The momentum of the recoil system used in decay dipole kinematics.
     * 25/05/2016 - Not currently used
     */
    //Lorentz5Momentum theRecoilMomentum;

    /**
     * The momentum of the recoil system after the splitting, used in decay dipole kinematics.
     */
    Lorentz5Momentum theSplitRecoilMomentum;



    int theOpenZBoundaries;

  protected:

    /**
     * Pointer to a check histogram object
     */
    Ptr<DipoleMCCheck>::ptr theMCCheck;

  private:

    /**
     * The static object used to initialize the description of this class.
     * Indicates that this is an abstract class.
     */
    static AbstractClassDescription<DipoleSplittingKinematics> initDipoleSplittingKinematics;

    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    DipoleSplittingKinematics & operator=(const DipoleSplittingKinematics &) = delete;

  };

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

  /** @cond TRAITSPECIALIZATIONS */

  /** This template specialization informs ThePEG about the
   *  base classes of DipoleSplittingKinematics. */
  template <>
  struct BaseClassTrait<Herwig::DipoleSplittingKinematics,1> {
    /** Typedef of the first base class of DipoleSplittingKinematics. */
    typedef HandlerBase NthBase;
  };

  /** This template specialization informs ThePEG about the name of
   *  the DipoleSplittingKinematics class and the shared object where it is defined. */
  template <>
  struct ClassTraits<Herwig::DipoleSplittingKinematics>
    : public ClassTraitsBase<Herwig::DipoleSplittingKinematics> {
    /** Return a platform-independent class name */
    static string className() { return "Herwig::DipoleSplittingKinematics"; }
    /**
     * The name of a file containing the dynamic library where the class
     * DipoleSplittingKinematics is implemented. It may also include several, space-separated,
     * libraries if the class DipoleSplittingKinematics depends on other classes (base classes
     * excepted). In this case the listed libraries will be dynamically
     * linked in the order they are specified.
     */
    static string library() { return "HwDipoleShower.so"; }
  };

  /** @endcond */

}

#endif /* HERWIG_DipoleSplittingKinematics_H */
