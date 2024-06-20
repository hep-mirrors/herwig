// -*- C++ -*-
#ifndef Herwig_DarkHadronSpectrum_H
#define Herwig_DarkHadronSpectrum_H
//
// This is the declaration of the DarkHadronSpectrum class.
//

#include "Herwig/Hadronization/HadronSpectrum.h"
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/PDT/StandardMatchers.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/PDT/EnumParticles.h>
#include "ThePEG/Repository/CurrentGenerator.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the DarkHadronSpectrum class.
 *
 * @see \ref DarkHadronSpectrumInterfaces "The interfaces"
 * defined for DarkHadronSpectrum.
 */
class DarkHadronSpectrum: public HadronSpectrum {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DarkHadronSpectrum(unsigned int opt);

  /**
   * The destructor.
   */
  virtual ~DarkHadronSpectrum();
  //@}

public:

  /** @name Partonic content */
  //@{

  /**
   * Return the id of the gluon
   */
  virtual long gluonId() const { return ParticleID::darkg; }

  /**
   * Return the ids of all hadronizing quarks
   */
  virtual const vector<long>& hadronizingQuarks() const {
    static vector<long> hadronizing = lightHadronizingQuarks();
    static vector<long> heavy = heavyHadronizingQuarks();
    hadronizing.insert(hadronizing.end(), heavy.begin(), heavy.end());
    return hadronizing;
  }

  /**
   * The light hadronizing quarks
   */
  virtual const vector<long>& lightHadronizingQuarks() const {
    if (long(_lightquarks.size()) != _nlightquarks) {
      for (long il=0; il<_nlightquarks; il++) {
        _lightquarks.push_back(il+_DarkHadOffset+1);
      }
    }
    return _lightquarks;
  }

  /**
   * The heavy hadronizing quarks
   */
  virtual const vector<long>& heavyHadronizingQuarks() const {
    if (long(_heavyquarks.size()) != _nheavyquarks) {
      for (long il=0; il<_nheavyquarks; il++) {
        _heavyquarks.push_back(il+_DarkHadOffset+1+_nlightquarks);
      }
    }
    return _heavyquarks;
  }

  /**
   * The lightest quarks, used for finding the lightest Hadron Pair
   */
  virtual const vector<long>& lightestQuarks() const {
    // May need to be updated in future for strange-like quarks
    return lightHadronizingQuarks();
  }


  /**
   * Return true if any of the possible three input particles contains
   * the indicated heavy quark.  false otherwise. In the case that
   * only the first particle is specified, it can be: an (anti-)quark,
   * an (anti-)diquark an (anti-)meson, an (anti-)baryon; in the other
   * cases, each pointer is assumed to be either (anti-)quark or
   * (anti-)diquark.
   */
  virtual bool hasHeavy(long, tcPDPtr, tcPDPtr = PDPtr(), tcPDPtr = PDPtr()) const {
    //ToDo: this should work for the heavyHadronizingQuarks
    return false;
  }

  //@}

  /**
   * Return the threshold for a cluster to split into a pair of hadrons.
   * This is normally the mass of the lightest hadron Pair, but can be
   * higher for heavy and exotic clusters
   */
  virtual Energy hadronPairThreshold(tcPDPtr par1, tcPDPtr par2) const;

  /**
   * Return the weight for the given flavour
   */
  virtual double pwtQuark(const long& id) const {
    return pwt(id);
  }

  /**
   * Dummy function as it will not be needed
   */
	virtual const vector<long>& lightHadronizingDiquarks() const {
		static vector<long> nothing;
		return nothing;
	};
  /**
   * The diquark weight.
   */
   double pwtDIquark() const {
    return _pwtDIquark;
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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   *
   *  The array _repwt is initialized using the interfaces to set different
   *  weights for different meson multiplets and the constructHadronTable()
   *  method called to complete the construction of the hadron tables.
   *
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

  
  /**
   * Return the id of the diquark (anti-diquark) made by the two 
   * quarks (antiquarks) of id specified in input (id1, id2).
   * Caller must ensure that id1 and id2 are quarks.
   */
  long makeDiquarkID(long id1, long id2, long pspin)  const;
  
  /**
   *  Weights for mesons
   */
  virtual double mesonWeight(long id) const;

  /**
   * Return true, if any of the possible input particle pointer is an exotic quark, e.g. Susy quark;
   * false otherwise.   
   */
  bool isExotic(tcPDPtr par1, tcPDPtr par2 = PDPtr(), tcPDPtr par3 = PDPtr())  const;

protected:

  /**
   *  Construct the table of hadron data
   *  This is the main method to initialize the hadron data (mainly the
   *  weights associated to each hadron, taking into account its spin, 
   *  eventual isoscalar-octect mixing, singlet-decuplet factor). This is
   *  the method that one should update when new or updated hadron data is
   *  available. 
   *
   *  This class implements the construction of the basic table but can be 
   *  overridden if needed in inheriting classes.
   *
   *  The rationale for factors used for diquarks involving different quarks can 
   *  be can be explained by taking a prototype example that in the  exact SU(2) limit,
   *  in which:
   *  \f[m_u=m_d\f] 
   *  \f[M_p=M_n=M_\Delta\f]
   *      and we will have equal numbers of u and d quarks produced.
   *      Suppose that we weight 1 the diquarks made of the same 
   *      quark and 1/2 those made of different quarks, the fractions
   *      of u and d baryons (p, n, Delta) we get are the following:
   *        - \f$\Delta^{++}\f$: 1 possibility only  u uu  with weight 1
   *        - \f$\Delta^-   \f$: 1 possibility only  d dd  with weight 1
   *        - \f$p,\Delta^+ \f$: 2 possibilities     u ud  with weight 1/2
   *                                                 d uu  with weight 1
   *        - \f$n,\Delta^0 \f$: 2 possibilities     d ud  with weight 1/2
   *                                                 u dd  with weight 1
   *      In the latter two cases, we have to take into account the 
   *      fact that  p  and  n  have spin 1/2 whereas  Delta+  and  Delta0
   *      have spin 3/2 therefore from phase space we get a double weight 
   *      for  Delta+  and  Delta0  relative to  p  and  n  respectively.
   *      Therefore the relative amount of these baryons that is
   *      produced is the following:
   *       # p = # n = ( 1/2 + 1 ) * 1/3 = 1/2
   *       # Delta++ = # Delta- = 1 = ( 1/2 + 1) * 2/3 # Delta+ = # Delta0
   *      which is correct, and therefore the weight 1/2 for the
   *      diquarks of different types of quarks is justified (at least
   *      in this limit of exact SU(2) ).
   */
  virtual void constructHadronTable();

  /**
   *  Access the parton weights
   */
   double pwt(long pid) const {
    map<long,double>::const_iterator it = _pwt.find(abs(pid));
    assert( it != _pwt.end() );
    return it->second;
  }

  /**
   *   Insert a meson in the table
   */
  virtual void insertMeson(HadronInfo a, int flav1, int flav2);

  /**
   * Methods for the mixing of \f$I=0\f$ mesons
   */
  //@{
  /**
   * Return the probability of mixing for Octet-Singlet isoscalar mixing,
   * the probability of the 
   * \f$\frac1{\sqrt{2}}(|u\bar{u}\rangle + |d\bar{d}\rangle)\f$ component
   * is returned.
   * @param angleMix The mixing angle in degrees (not radians)
   * @param order is 0 for no mixing, 1 for the first resonance of a pair,
   *                 2 for the second one.
   * The mixing is defined so that for example with \f$\eta-\eta'\f$ mixing where
   * the mixing angle is \f$\theta=-23^0$ with $\eta\f$ as the first particle
   * and \f$\eta'\f$ the second one.
   * The convention used is 
   * \f[\eta  = \cos\theta|\eta_{\rm octet  }\rangle
   *           -\sin\theta|\eta_{\rm singlet}\rangle\f]
   * \f[\eta' = \sin\theta|\eta_{\rm octet  }\rangle
   *           -\cos\theta|\eta_{\rm singlet}\rangle\f]
   * with 
   * \f[|\eta_{\rm singlet}\rangle = \frac1{\sqrt{3}}
   * \left[|u\bar{u}\rangle + |d\bar{d}\rangle +  |s\bar{s}\rangle\right]\f]
   * \f[|\eta_{\rm octet  }\rangle = \frac1{\sqrt{6}}
   * \left[|u\bar{u}\rangle + |d\bar{d}\rangle - 2|s\bar{s}\rangle\right]\f]
   */
   double probabilityMixing(const double angleMix,
				  const int order) const {
    static double convert=Constants::pi/180.0;
    if (order == 1)      
      return sqr( cos( angleMix*convert + atan( sqrt(2.0) ) ) );
    else if (order == 2) 
      return sqr( sin( angleMix*convert + atan( sqrt(2.0) ) ) );
    else                 
      return 1.;
  }

  /**
   * Returns the weight of given mixing state.
   * @param id The PDG code of the meson
   */
  virtual double mixingStateWeight(long id) const; 
  //@}

  virtual double specialQuarkWeight(double quarkWeight, long,
				    const Energy, tcPDPtr, tcPDPtr) const {
      return quarkWeight;
  }

  /**
   * The probability of producting a diquark.
   */
  double _pwtDIquark;

  /**
   * Singlet and Decuplet weights
   */
  //@{
  /**
   *  The singlet weight
   */
  double _sngWt; 

  /**
   *  The decuplet weight
   */
  double _decWt; 
  //@}

  /**
   * Return true if the two or three particles in input can be the components 
   * of a baryon; false otherwise.
   */
  virtual bool canBeBaryon(tcPDPtr par1, tcPDPtr par2 , tcPDPtr par3 = PDPtr())  const;

private:
  /**
   *  Option for the construction of the tables
   */ 
  unsigned int _topt;

  /**
   *  Which particles to produce for debugging purposes
   */
  unsigned int _trial;

  /**
   *  Prefix for Dark Hadron pdgID
   */
  int _DarkHadOffset = 4900000;

  /**
   *  The number of light quarks
   */
  int _nlightquarks;

  /**
   *  The number of heavy quarks
   */
  int _nheavyquarks;


  /**
   *  The pdgIds of the light quarks
   */
  mutable vector<long> _lightquarks = {};

  /**
   *  The pdgIds of the heavy quarks
   */
  mutable vector<long> _heavyquarks = {};

};

}

#endif /* Herwig_DarkHadronSpectrum_H */
