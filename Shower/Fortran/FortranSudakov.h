// -*- C++ -*-
#ifndef HERWIG_FortranSudakov_H
#define HERWIG_FortranSudakov_H
//
// This is the declaration of the FortranSudakov class.
//

#include "Herwig++/Shower/SplittingFunctions/SudakovFormFactor.h"
#include "Herwig++/Utilities/NewInterpolator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "FortranSudakov.fh"
#include "Herwig++/Shower/Couplings/ShowerAlphaQCD.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the FortranSudakov class.
 *
 * @see \ref FortranSudakovInterfaces "The interfaces"
 * defined for FortranSudakov.
 */
class FortranSudakov: public SudakovFormFactor {

public:

  /**
   * The default constructor.
   */
  inline FortranSudakov();

  /**
   *  Members to generate the scale of the next branching
   */
  //@{
  /**
   * Return the scale of the next time-like branching. If there is no 
   * branching then it returns Energy().
   * @param startingScale starting scale for the evolution
   * @param ids The PDG codes of the particles in the splitting
   * @param cc Whether this is the charge conjugate of the branching
   * defined.
   */
  virtual Energy generateNextTimeBranching(const Energy startingScale,
				           const IdList &ids,const bool cc);

  /**
   * Return the scale of the next space-like decay branching. If there is no 
   * branching then it returns Energy().
   * @param startingScale starting scale for the evolution
   * @param stoppingScale stopping scale for the evolution
   * @param minmass The minimum mass allowed for the spake-like particle.
   * @param ids The PDG codes of the particles in the splitting
   * @param cc Whether this is the charge conjugate of the branching
   * defined.
   */
  virtual Energy generateNextDecayBranching(const Energy startingScale,
					    const Energy stoppingScale,
					    const Energy minmass,
					    const IdList &ids,
					    const bool cc);

  /**
   * Return the scale of the next space-like branching. If there is no 
   * branching then it returns Energy().
   * @param startingScale starting scale for the evolution
   * @param ids The PDG codes of the particles in the splitting
   * @param x The fraction of the beam momentum
   * @param cc Whether this is the charge conjugate of the branching
   * defined.
   */
  virtual Energy generateNextSpaceBranching(const Energy startingScale,
		                            const IdList &ids,double x,
					    const bool cc);
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

  /**
   *  Return the shower cut-off for a particular particle
   */
  inline Energy cutOff(long) const;

  /**
   *  Find the tables for a given Sudakov
   * @param ids The particles
   * @param cc Whether or not these are hte conjugates of the particles defined
   */
  inline unsigned int findSudakov(const IdList & ids,const bool cc) const;

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
  static ClassDescription<FortranSudakov> initFortranSudakov;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FortranSudakov & operator=(const FortranSudakov &);

private:

  /**
   *  The order in \f$\alpha_S\f$ for the Sudakovs in the tables
   */
  unsigned int _sudord;

  /**
   *  The order of interpolation to use for the tables
   */
  unsigned int _inter;

  /**
   * Virtual mass cut-offs which are added to the constituent masses
   * in the parton shower
   */
  //@{
  /**
   * Quark cut-off
   */
  Energy _vqcut;

  /**
   * Gluon cut-off
   **/
  Energy _vgcut;

  /**
   * Photon cut-off
   */
  Energy _vpcut;
  //@}

  /**
   *  Number of enteries for the lookup tables
   */
  unsigned int _nqev;

  /**
   *  Vector of interpolators for the Sudakovs to interpolate
   *  to give the Sudakov given the scale
   */
  vector<NewInterpolatorPtr> _sudakovQ;

  /**
   *  Vector of interpolators for the Sudakovs to interpolate
   *  to give the scale given the Sudakov
   */
  vector<NewInterpolatorPtr> _sudakovP;

  /**
   *  Vector of interpolators for the Sudakovs to interpolate
   *  to give the Sudakov given the scale using linear extrapolation as a backup
   */
  vector<NewInterpolatorPtr> _linearQ;

  /**
   *  Vector of interpolators for the Sudakovs to interpolate
   *  to give the scale given the Sudakov using linear extrapolation as a backup
   */
  vector<NewInterpolatorPtr> _linearP;

};

}

namespace Herwig {

using namespace ThePEG;

/**
 *  This class provides the integrand for the Sudakov form factor
 */
class FortranSudakovIntegrand {

public:

  /**
   *  Constructor
   * @param qcdlam The value of \f$\Lambda_{\rm QCD}\f$ for three flavours
   * @param split  The splitting function to be integrated
   * @param ids    The particles involved in the splitting
   * @param sudord The order of \f$alpha_S\f$.
   * @param zord   Whether to use \f$z\f$ (false) or \f$1-z\f$ (true)
   *               in the splitting function
   * @param qmass Vector containing the constituent quark masses
   * @param alpha Pointer to the object calculating \f$\alpha_S\f$.
   */
  FortranSudakovIntegrand(Energy qcdlam, SplittingFnPtr split, IdList ids,
			  unsigned int sudord, bool zord, vector<Energy> qmass,
			  ShowerAlphaQCDPtr alpha);

  /**
   *  Return the value of the integrand
   */
  double operator() (double) const;

  /**
   *  Set the scales
   */
  inline void setScales(double qrat,double qlam);

protected:

  /**
   *  Integral of the 2nd order \f$alpha_S\f$ between two limits
   */
  inline double alphaIntegral(double lower, double upper,
			      unsigned int nf) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FortranSudakovIntegrand & operator=(const FortranSudakovIntegrand &);

  /**
   *   The value of \f$\Lambda_{\rm QCD}\f$ for three flavours
   */
  Energy _qcdlam;

  /**
   * The splitting function to be integrated
   */
  SplittingFnPtr _split; 

  /**
   *  The particles involved in the splitting
   */
  IdList _ids;

  /**
   * The order of \f$alpha_S\f$.
   */
  unsigned int _sudord;

  /**
   *  Whether to use \f$z\f$ (false) or \f$1-z\f$ (true) in the splitting function
   */
  bool _zord;

  /**
   *  Pointer to object calculating the strong coupling
   */
  ShowerAlphaQCDPtr _alpha;

  /**
   *  \f$b\f$ coefficients for the \f$\beta\f$ function.
   */
  vector<double> _bet;

  /**
   *  \f$b'\f$ coefficients for the \f$\beta\f$ function.
   */
  vector<double> _bep;

  /**
   *  Mass scales for the flavour thresholds
   */
  vector<Energy> _mumi,_muma;

  /**
   *  \f$\alpha_S\f$ values at the thresholds
   */
  vector<double> _almi,_alma;

  /**
   *  Integrals of \f$alpha_S\f$ between the limits
   */
  vector<double> _fint;
 
  /**
   *  Value of the scales
   */
  double _qrat,_qlam;
};
}







#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FortranSudakov. */
template <>
struct BaseClassTrait<Herwig::FortranSudakov,1> {
  /** Typedef of the first base class of FortranSudakov. */
  typedef Herwig::SudakovFormFactor NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FortranSudakov class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::FortranSudakov>
  : public ClassTraitsBase<Herwig::FortranSudakov> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::FortranSudakov"; }
  /**
   * The name of a file containing the dynamic library where the class
   * FortranSudakov is implemented. It may also include several, space-separated,
   * libraries if the class FortranSudakov depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "FortranSudakov.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "FortranSudakov.tcc"
#endif

#endif /* HERWIG_FortranSudakov_H */
