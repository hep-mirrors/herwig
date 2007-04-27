// -*- C++ -*-
#ifndef HERWIG_GeneralDipole_H
#define HERWIG_GeneralDipole_H
//
// This is the declaration of the GeneralDipole class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/CLHEPWrap/Lorentz5Vector.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "Herwig++/Utilities/Math.h"
#include "YFSFormFactors.h"
#include "GeneralDipole.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the GeneralDipole class.
 *
 * @see \ref GeneralDipoleInterfaces "The interfaces"
 * defined for GeneralDipole.
 */
class GeneralDipole: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline GeneralDipole();

  /**
   * The copy constructor.
   */
  inline GeneralDipole(const GeneralDipole &);

  /**
   * The destructor.
   */
  virtual ~GeneralDipole();
  //@}

public:

  /**
   * Member to generate the photons from the dipole
   * @param p The decaying particle
   * @param children The decay products
   * @return The decay products with additional radiation
   */
  virtual ParticleVector generatePhotons(const Particle & p,ParticleVector children);

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

protected
:
  /**
   * Return the photon multiplicity according to a Poissonian Distribution
   * with the supplied average
   * @param average The average
   * @return A value from tghe poisson distribution
   */
  inline int poisson(double average);

  /**
   * Average crude photon multiplicity for a given dipole
   * @param i First particle in the dipole
   * @param j Second particle in the dipole
   * @return The average photon multiplicity
   */
  inline double nbar(unsigned int i, unsigned int j);

  /**
   * Member which generates the photons
   */
  double makePhotons();

  /**
   * Generate the momenta of the photons 
   * @param i The first particle for the dipole
   * @param j The second particle for the dipole
   * @param nphoton The number of photons to produce
   */
  void photon(unsigned int i,unsigned int j,unsigned int nphoton);

  /**
   * Boost all the momenta from the dipole rest frame via the parent rest frame
   * to the lab
   * @param photons Whether or not there are any photons left to simplify the
   * calculations.
   * @return Whether or not it suceeded
   */
  bool boostMomenta(bool photons);

  /**
   *  Remove any photons which fail the energy cuts
   * @return Number of photons removed
   */
  unsigned int removePhotons();

  /**
   *  Reweights the dipole for the effects of the momentum rescaling
   */
  inline void reweightDipole();

  /**
   * The full YFS form factor
   * @param mom The momenta of the particles
   * @param mdipole The charge of the dipole
   * @param ecut The cut-off on the photon energy
   * @param poscharge Whether or not to include positively charged dipoles 
   */
  double YFSFormFactor(const vector<Lorentz5Momentum> & mom,
		       const vector<vector<Energy2> > mdipole,
		       const Energy & ecut,
		       bool poscharge) const;

  /**
   * Jacobian factor for the weight
   */
  inline double jacobianWeight();

  /**
   * Calculate the full weight for the dipoles
   */
  double fullDipoleWeight();

  /**
   * Weight for an individual dipole
   * @param ix First  particle in dipole
   * @param iy Second particle in dipole
   * @param iphot Photn
   * @return The weight
   */
  inline double dipoleWeight(unsigned int ix, unsigned int iy,unsigned int iphot);

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<GeneralDipole> initGeneralDipole;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GeneralDipole & operator=(const GeneralDipole &);

private:

  /**
   *   Momenta of the particles in the dipole rest frame
   */
  //@{
  /**
   *  Momenta of the charged particles in the dipole rest frame before radiation
   */
  vector<Lorentz5Momentum> _qdrf;

  /**
   *  Momenta of the charged particles in the dipole rest frame after radiation
   */
  vector<Lorentz5Momentum> _qnewdrf;

  /**
   *  Momenta of the photons in the dipole rest frame
   */
  vector<Lorentz5Momentum> _ldrf;

  /**
   * Total momentum of the photons in the dipole rest frame
   */
  Lorentz5Momentum _bigLdrf;
  //@}

  /**
   *  Momentum of the particles in the parent's rest frame
   */
  //@{
  /**
   *  Momenta of the charged particles in the parent's rest frame before radiation
   */
  vector<Lorentz5Momentum> _qprf;

  /**
   *  Momenta of the charged particles in the parent's rest frame after radiation
   */
  vector<Lorentz5Momentum> _qnewprf;

  /**
   *  Momenta of the photons in the parent rest frame
   */
  vector<Lorentz5Momentum> _lprf;

  /**
   * Total momentum of the photons in the parent rest frame
   */
  Lorentz5Momentum _bigLprf;
  //@}

  /**
   *  Momentum of the particles in the lab frame
   */
  //@{
  /**
   *  Momenta of the charged particles in the lab frame before radiation
   */
  vector<Lorentz5Momentum> _qlab;
  
  /**
   *  Momenta of the charged particles in the lab frame after  radiation
   */
  vector<Lorentz5Momentum> _qnewlab;

  /**
   *  Momenta of the photons in the lab frame
   */
  vector<Lorentz5Momentum> _llab;

  /**
   * Total momentum of the photons in the lab frame
   */
  Lorentz5Momentum _bigLlab;
  //@}

  /**
   *   Properties of the dipoles
   */
  //@{

  /**
   *  Masses of the different dipoles before radiation
   */
  vector<vector <Energy2> > _mdipole;

  /**
   *  Masses of the different dipoles after radiation
   */
  vector<vector <Energy2> > _mdipolenew;

  /**
   *  Opening angles for the different dipoles
   */
  vector<vector <double> > _cosij,_sinij;

  /**
   *  Charges of the different possible dipoles
   */
  vector< vector<int> > _zij;

  /**
   *  Average multiplicities for the different dipoles
   */
  vector<vector <double> > _nbar;
  //@}

  /**
   *  Properties of the decay products
   */
  //@{
  /**
   *  Boost from the rest frame to the lab frame
   */
  Hep3Vector _boosttolab;

  /**
   *  Number of decay products
   */
  unsigned int _nprod;

  /**
   *  Masses of the particles
   */
  vector<Energy> _m;

  /**
   *  Masses of the particles squared
   */
  vector<Energy2> _m2;

  /**
   *  Sum of the decay product masses
   */
  Energy _msum;

  /**
   *  The velocities of the particles before radiation
   */
  vector<double> _beta;

  /**
   * \f$1-\beta\f$ for the different particles before radiation
   */
  vector<double> _ombeta;

  /**
   *  The velocities of the particles after radiation
   */
  vector<double> _betanew;

  /**
   * \f$1-\beta\f$ for the different particles after radiation
   */
  vector<double> _ombetanew;
  //@}

  /**
   *  Limits on the photon energies
   */
  //@{
  /**
   *  The minimum photon energy in the boosted frame
   */
  Energy _emin;

  /**
   *  The minimum photon energy in the rest frame
   */
  Energy _eminrest;

  /**
   *  The minimum photon energy in the lab  frame
   */
  Energy _eminlab;

  /**
   *  The maximum photon energy
   */
  Energy _emax;
  //@}

  /**
   *   mass of the decay products after rescaling
   */
  Energy _roots;

  /**
   *  Rescaling factor for the momenta
   */
  double _rescalingfactor;

  /**
   *  Reweighting factors due to differences between the true and crude
   *  distributions
   */
  //@{
  /**
   *  Reweighting factor for the real emission
   */
  double _dipolewgt;

  /**
   *  Reweighting factor for the YFS form-factor
   */
  double _yfswgt;

  /**
   *  Reweighting factor due to phase space
   */
  double _jacobianwgt;

  /**
   *  Reweighting factor due to matrix element corrections
   */
  double _mewgt;

  /**
   *  Maximum weight
   */
  double _maxwgt;
  //@}

  /**
   *   Storage of the properties of the radiated photons
   */
  //@{
  /**
   *  Cosine of the photon angles
   */
  vector<double> _cosphot;

  /**
   *  Sine of the photon angles
   */
  vector<double> _sinphot;

  /**
   *  Full dipole weight for the photon before the momenta are adjusted
   */
  vector<double> _photonwgt;

  /**
   *  Whether a given photon passes the energy cut
   */
  vector<bool> _photcut;

  /**
   *  Emmiting particle for a given photon
   */
  vector<unsigned int> _photonemit;

  /**
   *  Spectator particle for a given photon
   */
  vector<unsigned int> _photonspect;
  //@}
  /**
   *  Maximum number of attempts to generate a result
   */
  unsigned int _maxtry;

  /**
   *  Photon multiplicity being generated
   */
  unsigned int _multiplicity;

  /**
   *  Parameters controlling the generation
   */
  //@{
  /**
   * Maximum number of photons to generate
   */
  unsigned int _nphotonmax;

  /**
   *  Type of unweighting to perform
   */
  unsigned int _mode;

  /**
   *  Option for the energy cut-off
   */
  unsigned int _energyopt;

  /**
   *  Option for the inclusion of higher order corrections
   */
  unsigned int _betaopt;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GeneralDipole. */
template <>
struct BaseClassTrait<Herwig::GeneralDipole,1> {
  /** Typedef of the first base class of GeneralDipole. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GeneralDipole class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::GeneralDipole>
  : public ClassTraitsBase<Herwig::GeneralDipole> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::GeneralDipole"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the GeneralDipole class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwDecRad.so"; }
};

/** @endcond */

}

#include "GeneralDipole.icc"

#endif /* HERWIG_GeneralDipole_H */
