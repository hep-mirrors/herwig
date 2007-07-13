// -*- C++ -*-
#ifndef HERWIG_HadronSelector_H
#define HERWIG_HadronSelector_H
//
// This is the declaration of the HadronSelector class.
//

#include "ThePEG/Interface/Interfaced.h"
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/PDT/StandardMatchers.h>
#include "HadronSelector.fh"

namespace Herwig {

using namespace ThePEG;
/**\ingroup Hadronization
 *  \class HadronSelector
 *  \brief This class selects the hadron flavours of a cluster decay.
 *  \author Philip Stephens
 *  \author Alberto Ribon
 *  \author Peter Richardson
 *
 *  This is the base class for the selection of either a pair of hadrons, or
 *  in some cases a single hadron. The different approaches which were previously
 *  implemented in this class are now implemented in the HwppSelector and Hw64Selector
 *  which inherit from this class.
 *
 *  This class implements a number of methods which are needed by all models
 *  and in addition contains the weights for the different meson multiplets and
 *  mixing of the light \f$I=0\f$ mesons.
 *
 * @see \ref HadronSelectorInterfaces "The interfaces"
 * defined for HadronSelector.
 * @see HwppSelector
 * @see Hw64Selector
 */
class HadronSelector: public Interfaced {

protected:

  /**
   *  The helper classes, definitions are now external
   */
  //@{
  /**
   *  The helper Kupco class
   */
  class Kupco;

  /**
   *  The helper HadronInfo class
   */
  class HadronInfo;

  /**
   * The type is used to contain all the hadrons info of a given flavour.
   */
  typedef set<HadronInfo> KupcoData;
  //@}

public:

  /**
   * The default constructor.
   */
  HadronSelector(unsigned int);

  /**
   * Method to return a pair of hadrons given the PDG codes of
   * two or three constituents
   * @param cluMass The mass of the cluster
   * @param id1 The PDG code of the first constituent
   * @param id2 The PDG code of the first constituent
   * @param id3 The PDG code of the first constituent
   */
  virtual pair<tcPDPtr,tcPDPtr> chooseHadronPair(const Energy cluMass, const long id1, 
						 const long id2, const long id3=0)
    throw(Veto, Stop, Exception) =0;

  /**
   * This returns the lightest pair of hadrons given by the flavours.
   *
   * Given the id of two (or three) constituents of a cluster, it returns
   * the id of the two lightest hadrons with proper flavour numbers.
   * Furthermore, the first of the two hadrons must have the constituent with
   * id1, and the second must have the constituent with id2. 
   * At the moment it does *nothing* in the case that also id3 is present.
   *
   * The method is implemented by calling twice lightestHadron, 
   * once with (id1.-idPartner) , and once with (id2,idPartner)
   * where  idPartner (with sign defined properly) is either the id of
   * d or u . In fact, the idea is that whatever the flavour of id1 
   * and id2, no matter if (anti-)quark or (anti-)diquark, the lightest
   * pair of hadrons containing flavour id1 and id2 will have either 
   * flavour d or u, being the lightest quarks.
   * The method returns the pair (0,0) if anything goes wrong. 
   *
   * The method assumes id3==0 (otherwise we don't know how to proceed: a 
   * possible, trivial way would be to randomly select two of the three 
   * (anti-)quarks and treat them as a (anti-)diquark, reducing the problem
   * to two components as treated below.
   * In the normal (two components) situation, the strategy is the following:
   * treat in the same way the two possibilities:  (d dbar)  (i=0) and  
   * (u ubar)  (i=1)  as the pair quark-antiquark necessary to form a
   * pair of hadrons containing the input flavour  id1  and  id2; finally,
   * select the one that produces the lightest pair of hadrons, compatible
   * with the charge conservation contraint.
   */
  pair<tcPDPtr,tcPDPtr> lightestHadronPair(const long id1, const long id2,
					   const long id3=0) const;

  /**
   *  Returns the mass of the lightest pair of hadrons with the given ids.
   * @param id1 The PDG code of the first constituent
   * @param id2 The PDG code of the first constituent
   * @param id3 The PDG code of the first constituent
   */
  inline Energy massLightestHadronPair(const long id1, const long id2,
				       const long id3=0) const;

  /**
   * Returns the lightest hadron formed by the given ids.
   *
   * Given the id of two (or three) constituents of a cluster, it returns
   * the id of the lightest hadron with proper flavour numbers.
   * At the moment it does *nothing* in the case that also id3 is present.
   * @param id1 The PDG code of the first constituent
   * @param id2 The PDG code of the first constituent
   * @param id3 The PDG code of the first constituent
   */
  inline tcPDPtr lightestHadron(const long id1, const long id2,
				const long id3=0) const;

  /**
   * Return the nominal mass of the hadron with id returned by lightestHadron()
   * @param id1 The PDG code of the first constituent
   * @param id2 The PDG code of the first constituent
   * @param id3 The PDG code of the first constituent
   */
  inline Energy massLightestHadron(const long id1, const long id2,
				   const long id3=0) const;

  /**
   *  Returns the mass of the lightest pair of baryons with the given ids.
   * @param id1 The PDG code of the first constituent
   * @param id2 The PDG code of the first constituent
   */
  inline Energy massLightestBaryonPair(const long id1, const long id2) const;

  /**
   *  Return the weights for the different quarks and diquarks
   */
  //@{
  /**
   * The down quark weight.
   */
  inline double pwtDquark()  const; 

  /**
   * The up quark weight.
   */
  inline double pwtUquark()  const; 

  /**
   * The strange quark weight.
   */
  inline double pwtSquark()  const; 

  /**
   * The charm quark weight.
   */
  inline double pwtCquark()  const; 

  /**
   * The bottom quark weight.
   */
  inline double pwtBquark()  const; 

  /**
   * The diquark weight.
   */
  inline double pwtDIquark() const; 
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
  virtual void doinit() throw(InitException);
  //@}

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
   *  Access to the table of hadrons
   */
  inline map<pair<long,long>,KupcoData> & table();

  /**
   *  Access to the list of partons
   */
  inline vector<long> & partons();

  /**
   *  Access the parton weights
   */
  inline map<long,double> & pwt();

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
  inline double probabilityMixing(const double angleMix, const int order);

  /**
   * Returns the weight of given mixing state.
   * @param id The PDG code of the meson
   */
  virtual double mixingStateWeight(long id); 
  //@}

  /**
   * Calculates a special weight specific to  a given hadron.
   * @param id The PDG code of the hadron
   */
  double specialWeight(long id);

  /**
   * This method returns the proper sign ( > 0 hadron; < 0 anti-hadron )
   * for the input PDG id  idHad > 0, suppose to be made by the
   * two constituents of ids: idQ1 and idQ2 (both with proper sign).
   * In the case of failure, it returns 0.
   */
  int  signHadron(const int idQ1, const int idQ2, const tcPDPtr hadron) const;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<HadronSelector> initHadronSelector;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HadronSelector & operator=(const HadronSelector &);

private:

  /**
   *  The PDG codes of the constituent particles allowed
   */
  vector<long> _partons;

  /**
   *  Map of the PDG codes of the partons and their charges to avoid
   *  using getParticleData too much
   */
  map<long,long> _charge;

  /**
   *  The PDG codes of the hadrons which cannot be produced in the hadronization
   */
  vector<long> _forbidden;

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
   *  The weights for the excited meson multiplets
   */
  vector<vector<vector<double> > > _repwt;

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
   * The table of hadron data
   */
  map<pair<long,long>,KupcoData> _table;

  /**
   * Enums so arrays can be statically allocated
   */
  //@{
  /**
   * Defines values for array sizes. L,J,N max values for excited mesons.
   */
  enum MesonMultiplets { Lmax = 3, Jmax = 4, Nmax = 4}; 
  //@}

  /**
   *  Option for the construction of the tables
   */ 
  unsigned int _topt;

  /**
   *  Which particles to produce for debugging purposes
   */
  unsigned int _trial;
};

/** \ingroup Hadronization
 *  \class HadronInfo 
 *  \brief Class used to store all the hadron information for easy access.
 *  \author Philip Stephens
 *  
 *  Note that:
 *  - the hadrons in _table can be filled in any ordered 
 *    w.r.t. the mass value, and flavours for different
 *    groups (for instance, (u,s) hadrons don't need to
 *    be placed after (d,s) or any other flavour), but 
 *    all hadrons with the same flavours must be consecutive
 *    ( for instance you cannot alternate hadrons of type
 *    (d,s) with those of flavour (u,s) ).
 *    Furthermore, it is assumed that particle and antiparticle
 *    have the same weights, and therefore only one of them
 *    must be entered in the table: we have chosen to refer
 *    to the particle, defined as PDG id > 0, although if
 *    an anti-particle is provided in input it is automatically
 *    transform to its particle, simply by taking the modulus
 *    of its id.
 */
class Herwig::HadronSelector::HadronInfo {  

public:

  /**
   *  Constructor
   * @param idin The PDG code of the hadron
   * @param datain The pointer to the ParticleData object
   * @param swtin  The singlet/decuplet/orbital factor
   * @param massin The mass of the hadron
   */
  HadronInfo(long idin=0,tPDPtr datain=tPDPtr(),double swtin=1.,Energy massin=0.*MeV);

  /**
   *  Comparision operator on mass
   */
  inline bool operator<(const HadronInfo &x) const;

  /**
   * The hadrons id.
   */
  long  id;
  
  /**
   * pointer to ParticleData, to get the spin, etc...
   */
  tPDPtr ptrData;
  
  /**
   * singlet/decuplet/orbital factor 
   */
  double swtef;
  
  /**
   * mixing factor
   */
  double wt;          
  
  /**
   * (2*J+1)*wt*swtef
   */
  double overallWeight;
  
  /**
   * The hadrons mass
   */
  Energy mass;

  /**
   *  Rescale the weight for a given hadron
   */
  void rescale(double x) const { 
    const_cast<HadronInfo*>(this)->overallWeight *= x; 
  }

  /**
   * Friend method used to print the value of a table element.
   */
  friend PersistentOStream & operator<< (PersistentOStream & os, 
					 const HadronInfo & hi ) {
    os << hi.id << hi.ptrData << hi.swtef << hi.wt << hi.overallWeight << ounit(hi.mass,GeV);
    return os;
  }
  
  /**
   * Friend method used to read in the value of a table element.
   */
  friend PersistentIStream & operator>> (PersistentIStream & is, 
					 HadronInfo & hi ) {
    is >> hi.id >> hi.ptrData >> hi.swtef >> hi.wt >> hi.overallWeight >> iunit(hi.mass,GeV);
    return is;
  }
};

/** \ingroup Hadronization
 *  \class Kupco
 *  \brief Class designed to make STL routines easy to use.
 *  \author Philip Stephens
 *
 *  This class is used to generate a list of the hadron pairs which can 
 *  be produced that allows easy traversal and quick access.
 */
class Herwig::HadronSelector::Kupco {

public:
  
  /**
   *  Constructor
   * @param inidQ PDG code of the quark drawn from the vacuum.
   * @param inhad1 ParticleData for the first hadron produced.
   * @param inhad2 ParticleData for the second hadron produced.
   * @param inwgt  The weight for the hadron pair 
   */
  inline Kupco(long inidQ,tcPDPtr inhad1,tcPDPtr inhad2, Energy inwgt);
  
  /**
   * id of the quark drawn from the vacuum.
   */
  long idQ;
   
  /**
   * The ParticleData object for the first hadron produced.
   */
  tcPDPtr hadron1;
  
  /**
   * The ParticleData object for the second hadron produced.
   */
  tcPDPtr hadron2;

  /**
   * Weight factor of this componation.
   */
  Energy weight;
};
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HadronSelector. */
template <>
struct BaseClassTrait<Herwig::HadronSelector,1> {
  /** Typedef of the first base class of HadronSelector. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HadronSelector class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HadronSelector>
  : public ClassTraitsBase<Herwig::HadronSelector> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::HadronSelector"; }
};

/** @endcond */

}

#include "HadronSelector.icc"

#endif /* HERWIG_HadronSelector_H */
