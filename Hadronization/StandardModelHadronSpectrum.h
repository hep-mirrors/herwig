// -*- C++ -*-
#ifndef Herwig_StandardModelHadronSpectrum_H
#define Herwig_StandardModelHadronSpectrum_H
//
// This is the declaration of the StandardModelHadronSpectrum class.
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
 * Here is the documentation of the StandardModelHadronSpectrum class.
 *
 * @see \ref StandardModelHadronSpectrumInterfaces "The interfaces"
 * defined for StandardModelHadronSpectrum.
 */
class StandardModelHadronSpectrum: public HadronSpectrum {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  StandardModelHadronSpectrum(unsigned int opt);

  /**
   * The destructor.
   */
  virtual ~StandardModelHadronSpectrum();
  //@}

public:

  /** @name Partonic content */
  //@{

  /**
   * Return the id of the gluon
   */
  virtual long gluonId() const { return ParticleID::g; }

  /**
   * Return the ids of all hadronizing quarks
   */
  virtual const vector<long>& hadronizingQuarks() const {
    static vector<long> hadronizing =
      { ParticleID::d, ParticleID::u, ParticleID::s, ParticleID::c, ParticleID::b };
    return hadronizing;
  }

  /**
   * The light hadronizing quarks
   */
  virtual const vector<long>& lightHadronizingQuarks() const {
    static vector<long> light =
      { ParticleID::d, ParticleID::u, ParticleID::s };
    return light;
  }

  /**
   * The heavy hadronizing quarks
   */
  virtual const vector<long>& heavyHadronizingQuarks() const {
    static vector<long> heavy =
      { ParticleID::c, ParticleID::b };
    return heavy;
  }

  /**
   * The lightest quarks, used for finding the lightest Hadron Pair
   */
  virtual const vector<long>& lightestQuarks() const {
    static vector<long> light =
      { ParticleID::d, ParticleID::u};
    return light;
  }

  /**
   * Return true if any of the possible three input particles contains
   * the indicated heavy quark.  false otherwise. In the case that
   * only the first particle is specified, it can be: an (anti-)quark,
   * an (anti-)diquark an (anti-)meson, an (anti-)baryon; in the other
   * cases, each pointer is assumed to be either (anti-)quark or
   * (anti-)diquark.
   */
  virtual bool hasHeavy(long id, tcPDPtr par1, tcPDPtr par2 = PDPtr(), tcPDPtr par3 = PDPtr()) const {
    if ( abs(id) == ParticleID::c )
      return hasCharm(par1,par2,par3);
    if ( abs(id) == ParticleID::b )
      return hasBottom(par1,par2,par3);
    return false;
  }

  //@}

  /**
   * Return the quark flavour which should be considered to set the
   * minimum mass for a minimal cluster splitting.
   */
  virtual long minimalSplitQuark() const {
    //cout << "ParticleID::d= " << ParticleID::d << endl;
    //cout << getParticleData(ParticleID::d)->id() << endl;
    return ParticleID::d;
  }

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
    switch(id) {
    case ParticleID::d: return pwtDquark(); break;
    case ParticleID::u: return pwtUquark(); break;
    case ParticleID::s: return pwtSquark(); break;
    case ParticleID::c: return pwtCquark(); break;
    case ParticleID::b: return pwtBquark(); break;
    }
    return 0.;
  }

  /**
   * The down quark weight.
   */
   double pwtDquark()  const {
    return _pwtDquark;
  } 

  /**
   * The up quark weight.
   */
   double pwtUquark()  const { 
    return _pwtUquark;
  }

  /**
   * The strange quark weight.
   */
   double pwtSquark()  const { 
    return _pwtSquark;
  }

  /**
   * The charm quark weight.
   */
   double pwtCquark()  const { 
    return _pwtCquark;
  }

  /**
   * The bottom quark weight.
   */
   double pwtBquark()  const { 
    return _pwtBquark;
  } 
  
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
  long makeDiquarkID(long id1, long id2)  const;

  /**
   * Return true if any of the possible three input particles has
   * b-flavour; 
   * false otherwise. In the case that only the first particle is specified,
   * it can be: an (anti-)quark, an (anti-)diquark
   * an (anti-)meson, an (anti-)baryon; in the other cases, each pointer
   * is assumed to be either (anti-)quark or (anti-)diquark.
   */
  bool hasBottom(tcPDPtr par1, tcPDPtr par2 = PDPtr(), tcPDPtr par3 = PDPtr())  const;
  /**
   * Return true if any of the possible three input particles has 
   * c-flavour; 
   * false otherwise.In the case that only the first pointer is specified,
   * it can be: a (anti-)quark, a (anti-)diquark
   * a (anti-)meson, a (anti-)baryon; in the other cases, each pointer
   * is assumed to be either (anti-)quark or (anti-)diquark.
   */
  bool hasCharm(tcPDPtr par1, tcPDPtr par2 = PDPtr(), tcPDPtr par3 = PDPtr())  const;
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
  
private:

  /**
   *  The weights for the different quarks and diquarks
   */
  //@{
  /**
   * The probability of producting a down quark.
   */
  double _pwtDquark;

  /**
   * The probability of producting an up quark.
   */
  double _pwtUquark;

  /**
   * The probability of producting a strange quark.
   */
  double _pwtSquark;

  /**
   * The probability of producting a charm quark.
   */
  double _pwtCquark;

  /**
   * The probability of producting a bottom quark.
   */
  double _pwtBquark;

  /**
   * The probability of producting a diquark.
   */
  double _pwtDIquark;

  /**
   * Weights for quarks and diquarks.
   */
  map<long,double> _pwt;
  //@}

  /**
   *  The mixing angles for the \f$I=0\f$ mesons containing light quarks
   */
  //@{
  /**
   *  The \f$\eta-\eta'\f$ mixing angle 
   */
  double _etamix;

  /**
   *  The \f$\phi-\omega\f$ mixing angle
   */
  double _phimix;

  /**
   *  The \f$h_1'-h_1\f$ mixing angle
   */
  double _h1mix;

  /**
   *  The \f$f_0(1710)-f_0(1370)\f$ mixing angle
   */
  double _f0mix;

  /**
   *  The \f$f_1(1420)-f_1(1285)\f$ mixing angle
   */
  double _f1mix;

  /**
   *  The \f$f'_2-f_2\f$ mixing angle
   */
  double _f2mix;

  /**
   *  The \f$\eta_2(1870)-\eta_2(1645)\f$ mixing angle
   */
  double _eta2mix;

  /**
   *  The \f$\phi(???)-\omega(1650)\f$ mixing angle
   */
  double _omhmix;

  /**
   *  The \f$\phi_3-\omega_3\f$ mixing angle
   */
  double _ph3mix;

  /**
   *  The \f$\eta(1475)-\eta(1295)\f$ mixing angle
   */
  double _eta2Smix;

  /**
   *  The \f$\phi(1680)-\omega(1420)\f$ mixing angle
   */
  double _phi2Smix;
  //@}

  /**
   *  The weights for the various meson multiplets to be used to supress the
   * production of particular states
   */
  //@{
  /**
   *  The weights for the \f$\phantom{1}^1S_0\f$ multiplets
   */
  vector<double> _weight1S0;

  /**
   *  The weights for the \f$\phantom{1}^3S_1\f$ multiplets
   */
  vector<double> _weight3S1;

  /**
   *  The weights for the \f$\phantom{1}^1P_1\f$ multiplets
   */
  vector<double> _weight1P1;

  /**
   *  The weights for the \f$\phantom{1}^3P_0\f$ multiplets
   */
  vector<double> _weight3P0;

  /**
   *  The weights for the \f$\phantom{1}^3P_1\f$ multiplets
   */
  vector<double> _weight3P1;

  /**
   *  The weights for the \f$\phantom{1}^3P_2\f$ multiplets
   */
  vector<double> _weight3P2;

  /**
   *  The weights for the \f$\phantom{1}^1D_2\f$ multiplets
   */
  vector<double> _weight1D2;

  /**
   *  The weights for the \f$\phantom{1}^3D_1\f$ multiplets
   */
  vector<double> _weight3D1;

  /**
   *  The weights for the \f$\phantom{1}^3D_2\f$ multiplets
   */
  vector<double> _weight3D2;

  /**
   *  The weights for the \f$\phantom{1}^3D_3\f$ multiplets
   */
  vector<double> _weight3D3;
  //@}

  /**
   *  Option for the construction of the tables
   */ 
  unsigned int _topt;

  /**
   *  Which particles to produce for debugging purposes
   */
  unsigned int _trial;

  /**
   * @name A parameter used for determining when clusters are too light.
   *
   * This parameter is used for setting the lower threshold, \f$ t \f$ as
   * \f[ t' = t(1 + r B^1_{\rm lim}) \f]
   * where \f$ r \f$ is a random number [0,1].
   */
  //@{
  double _limBottom;
  double _limCharm;
  double _limExotic;
  //@}

};

}

#endif /* Herwig_StandardModelHadronSpectrum_H */
