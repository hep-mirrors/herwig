// -*- C++ -*-
//
// HadronSelector.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
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
#include <ThePEG/Repository/EventGenerator.h>
#include "HadronSelector.fh"
#include "HadronInfo.h"
#include "Kupco.h"

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
 *  in some cases a single hadron. The different approaches which were
 *  previously implemented in this class are now implemented in the
 *  HwppSelector and Hw64Selector which inherit from this class.
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

public:

  /**
   * The default constructor.
   */
  HadronSelector(unsigned int);

  /**
   * Method to return a pair of hadrons given the PDG codes of
   * two or three constituents
   * @param cluMass The mass of the cluster
   * @param par1 The first constituent
   * @param par2 The second constituent
   * @param par3 The third constituent
   */
  virtual pair<tcPDPtr,tcPDPtr> chooseHadronPair(const Energy cluMass, 
						 tcPDPtr par1, tcPDPtr par2) const;

  /**
   * Select the single hadron for a cluster decay
   * return null pointer if not a single hadron decay
   * @param par1 1st constituent
   * @param par2 2nd constituent
   * @param mass Mass of the cluster
   */
  tcPDPtr chooseSingleHadron(tcPDPtr par1, tcPDPtr par2, Energy mass) const;

  /**
   * This returns the lightest pair of hadrons given by the flavours.
   *
   * Given the two (or three) constituents of a cluster, it returns
   * the two lightest hadrons with proper flavour numbers.
   * Furthermore, the first of the two hadrons must have the constituent with
   * par1, and the second must have the constituent with par2.
   * \todo At the moment it does *nothing* in the case that also par3 is present.
   *
   * The method is implemented by calling twice lightestHadron,
   * once with (par1,quarktopick->CC()) ,and once with (par2,quarktopick)
   * where quarktopick is either the pointer to
   * d or u quarks . In fact, the idea is that whatever the flavour of par1
   * and par2, no matter if (anti-)quark or (anti-)diquark, the lightest
   * pair of hadrons containing flavour par1 and par2 will have either
   * flavour d or u, being the lightest quarks.
   * The method returns the pair (PDPtr(),PDPtr()) if anything goes wrong.
   *
   * \todo The method assumes par3 == PDPtr() (otherwise we don't know how to proceed: a
   * possible, trivial way would be to randomly select two of the three
   * (anti-)quarks and treat them as a (anti-)diquark, reducing the problem
   * to two components as treated below.
   * In the normal (two components) situation, the strategy is the following:
   * treat in the same way the two possibilities:  (d dbar)  (i=0) and
   * (u ubar)  (i=1)  as the pair quark-antiquark necessary to form a
   * pair of hadrons containing the input flavour  par1  and  par2; finally,
   * select the one that produces the lightest pair of hadrons, compatible
   * with the charge conservation constraint.
   */
  pair<tcPDPtr,tcPDPtr> lightestHadronPair(tcPDPtr ptr1, tcPDPtr ptr2,
					   tcPDPtr ptr3 = PDPtr ()) const;

  /**
   *  Returns the mass of the lightest pair of hadrons with the given particles
   * @param ptr1 is the first  constituent
   * @param ptr2 is the second constituent
   * @param ptr3 is the third  constituent
   */
    Energy massLightestHadronPair(tcPDPtr ptr1, tcPDPtr ptr2,
				  tcPDPtr ptr3 = PDPtr ()) const  {
    pair<tcPDPtr,tcPDPtr> pairData = lightestHadronPair(ptr1, ptr2, ptr3);
    if ( ! pairData.first || ! pairData.second ) return ZERO;
    return ( pairData.first->mass() + pairData.second->mass() );
  }

  /**
   * Returns the lightest hadron formed by the given particles.
   *
   * Given the id of two (or three) constituents of a cluster, it returns
   * the  lightest hadron with proper flavour numbers.
   * At the moment it does *nothing* in the case that also 'ptr3' present.
   * @param ptr1 is the first  constituent
   * @param ptr2 is the second constituent
   * @param ptr3 is the third  constituent
   */
   tcPDPtr lightestHadron(tcPDPtr ptr1, tcPDPtr ptr2,
			  tcPDPtr ptr3 = PDPtr ()) const;

  /**
   * Returns the hadrons below the constituent mass threshold formed by the given particles,
   * together with their total weight
   *
   * Given the id of two (or three) constituents of a cluster, it returns
   * the  lightest hadron with proper flavour numbers.
   * At the moment it does *nothing* in the case that also 'ptr3' present.
   * @param threshold The theshold
   * @param ptr1 is the first  constituent
   * @param ptr2 is the second constituent
   * @param ptr3 is the third  constituent
   */
  vector<pair<tcPDPtr,double> >
  hadronsBelowThreshold(Energy threshold,
			tcPDPtr ptr1, tcPDPtr ptr2,
			tcPDPtr ptr3 = PDPtr ()) const;

  /**
   * Return the nominal mass of the hadron returned by lightestHadron()
   * @param ptr1 is the first  constituent
   * @param ptr2 is the second constituent
   * @param ptr3 is the third  constituent
   */
   Energy massLightestHadron(tcPDPtr ptr1, tcPDPtr ptr2,
#ifndef NDEBUG
  				   tcPDPtr ptr3 = PDPtr ()) const {
#else
                                   tcPDPtr = PDPtr ()) const {
#endif
    // The method assumes ptr3 == empty
    assert(!(ptr3));
    // find entry in the table
    pair<long,long> ids(abs(ptr1->id()),abs(ptr2->id()));
    HadronTable::const_iterator tit=_table.find(ids);
    // throw exception if flavours wrong
    if(tit==_table.end()||tit->second.empty())
      throw Exception() <<  "HadronSelector::massLightestHadron "
			<< "failed for particle" << ptr1->id()  << " "
			<< ptr2->id()
			<< Exception::eventerror;
    // return the mass
    return tit->second.begin()->mass;
  }

  /**
   *  Returns the mass of the lightest pair of baryons.
   * @param ptr1 is the first  constituent
   * @param ptr2 is the second constituent
   */
  Energy massLightestBaryonPair(tcPDPtr ptr1, tcPDPtr ptr2) const;

  /**
   *  Access the parton weights
   */
   double pwt(long pid) const {
    map<long,double>::const_iterator it = _pwt.find(abs(pid));
    assert( it != _pwt.end() );
    return it->second;
  }

  /**
   *  Force baryon/meson selection
   */
  virtual pair<bool,bool> selectBaryon(const Energy cluMass, tcPDPtr par1, tcPDPtr par2) const;

  /**
   *  Strange quark weight
   */
  virtual double strangeWeight(const Energy cluMass, tcPDPtr par1, tcPDPtr par2) const;

  /**
   * Return the particle data of the diquark (anti-diquark) made by the two
   * quarks (antiquarks) par1, par2.
   * @param par1 (anti-)quark data pointer
   * @param par2 (anti-)quark data pointer
   */
  virtual PDPtr makeDiquark(tcPDPtr par1, tcPDPtr par2);
  
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

protected:

  /**
   *  A sub-function of HadronSelector::constructHadronTable().
   *  It receives the information of a prospective Hadron and inserts it
   *  into the hadron table construct.
   *  @param particle is a particle data pointer to the hadron
   *  @param flav1 is the first  constituent of the hadron
   *  @param flav2 is the second constituent of the hadron
   */
  void insertToHadronTable(tPDPtr &particle, int flav1, int flav2);

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
  const HadronTable & table() const {
    return _table;
  }

  /**
   *  Access to the table of hadrons
   */
  HadronTable & table() {
    return _table;
  }
  
  /**
   *  Access to the list of partons
   */
  const vector<PDPtr> & partons() const {
    return _partons;
  }
  
  /**
   *  Access to the list of partons
   */
  vector<PDPtr> & partons() {
    return _partons;
  }

  /**
   *  Access the parton weights
   */
  map<long,double> & pwt() {
    return _pwt;
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

  /**
   * Calculates a special weight specific to  a given hadron.
   * @param id The PDG code of the hadron
   */
  double specialWeight(long id) const {
    const int pspin = id % 10;
    // Only K0L and K0S have pspin == 0, should
    // not get them until Decay step
    assert( pspin != 0 );
    // Baryon : J = 1/2 or 3/2
    if(pspin%2==0)
      return baryonWeight(id);
    // Meson
    else 
      return mesonWeight(id); 
  }
  
  /**
   *  Weights for mesons
   */
  virtual double mesonWeight(long id) const;

  /**
   *  Weights for baryons
   */
  virtual double baryonWeight(long id) const = 0;

  /**
   * This method returns the proper sign ( > 0 hadron; < 0 anti-hadron )
   * for the input PDG id  idHad > 0, suppose to be made by the
   * two constituent particle pointers: par1 and par2 (both with proper sign).
   */
  int signHadron(tcPDPtr ptr1, tcPDPtr ptr2, tcPDPtr hadron) const;

  /**
   *   Insert a meson in the table
   */
  virtual void insertMeson(HadronInfo a, int flav1, int flav2);

  /**
   *   Insert a spin\f$\frac12\f$ baryon in the table
   */
  virtual void insertOneHalf(HadronInfo a, int flav1, int flav2);

  /**
   *   Insert a spin\f$\frac32\f$ baryon in the table
   */
  virtual void insertThreeHalf(HadronInfo a, int flav1, int flav2);
  
private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HadronSelector & operator=(const HadronSelector &) = delete;

private:

  /**
   *  The PDG codes of the constituent particles allowed
   */
  vector<PDPtr> _partons;

  /**
   *  The PDG codes of the hadrons which cannot be produced in the hadronization
   */
  vector<PDPtr> _forbidden;

  /**
   * Weights for quarks and diquarks.
   */
  map<long,double> _pwt;

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
   * The table of hadron data
   */
  HadronTable _table;

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

  /**
   *  Option for the selection of hadrons below the pair threshold
   */
  unsigned int belowThreshold_;
};


}

#endif /* HERWIG_HadronSelector_H */
