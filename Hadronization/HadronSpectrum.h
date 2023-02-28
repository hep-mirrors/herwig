// -*- C++ -*-
#ifndef Herwig_HadronSpectrum_H
#define Herwig_HadronSpectrum_H
//
// This is the declaration of the HadronSpectrum class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "HadronSpectrum.fh"
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/ParticleData.h>
#include "Kupco.h"
/* These last two imports don't seem to be used here, but are needed for other
    classes which import this. Should tidy up at some point*/
#include <ThePEG/PDT/StandardMatchers.h>
#include <ThePEG/Repository/EventGenerator.h>

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the HadronSpectrum class.
 *
 * @see \ref HadronSpectrumInterfaces "The interfaces"
 * defined for HadronSpectrum.
 */
class HadronSpectrum: public Interfaced {

public:

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
  class HadronInfo {  

  public:
    
    /**
     *  Constructor
     * @param idin The PDG code of the hadron
     * @param datain The pointer to the ParticleData object
     * @param swtin  The singlet/decuplet/orbital factor
     * @param massin The mass of the hadron
     */
    HadronInfo(long idin=0, tPDPtr datain=tPDPtr(),
	       double swtin=1., Energy massin=ZERO)
      : id(idin), ptrData(datain), swtef(swtin), wt(1.0), overallWeight(0.0),
	mass(massin)
    {}
    
    /**
     *  Comparision operator on mass
     */
     bool operator<(const HadronInfo &x) const {
      if(mass!=x.mass) return mass < x.mass;
      else return id < x.id;
    }
    
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
      os << hi.id << hi.ptrData << hi.swtef << hi.wt 
	 << hi.overallWeight << ounit(hi.mass,GeV);
      return os;
    }
    
    /**
     * debug output
     */
    friend ostream & operator<< (ostream & os, const HadronInfo & hi ) {
      os << std::scientific << std::showpoint
	 << std::setprecision(4)
	 << setw(2)
	 << hi.id << '\t'
	 << hi.swtef << '\t'
	 << hi.wt << '\t'
	 << hi.overallWeight << '\t'
	 << ounit(hi.mass,GeV);
      return os;
    }
    
    /**
     * Friend method used to read in the value of a table element.
     */
    friend PersistentIStream & operator>> (PersistentIStream & is, 
					   HadronInfo & hi ) {
      is >> hi.id >> hi.ptrData >> hi.swtef >> hi.wt 
	 >> hi.overallWeight >> iunit(hi.mass,GeV);
      return is;
    }
  };
  

public:

  /**
   *  The helper classes
   */
  //@{
  /**
   * The type is used to contain all the hadrons info of a given flavour.
   */
  typedef set<HadronInfo> KupcoData;
  //@}

  /**
   * The hadron table type.
   */
  typedef map<pair<long,long>,KupcoData> HadronTable;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  HadronSpectrum();

  /**
   * The destructor.
   */
  virtual ~HadronSpectrum();
  //@}

public:

  /** @name Partonic content */
  //@{

  /**
   * Return the id of the gluon
   */
  virtual long gluonId() const = 0;

  /**
   * Return the ids of all hadronizing quarks
   */
  virtual const vector<long>& hadronizingQuarks() const = 0;

  /**
   * The light hadronizing quarks
   */
  virtual const vector<long>& lightHadronizingQuarks() const = 0;

  /**
   * The heavy hadronizing quarks
   */
  virtual const vector<long>& heavyHadronizingQuarks() const = 0;

  /**
   * Return true if any of the possible three input particles contains
   * the indicated heavy quark.  false otherwise. In the case that
   * only the first particle is specified, it can be: an (anti-)quark,
   * an (anti-)diquark an (anti-)meson, an (anti-)baryon; in the other
   * cases, each pointer is assumed to be either (anti-)quark or
   * (anti-)diquark.
   */
  virtual bool hasHeavy(long id, tcPDPtr par1, tcPDPtr par2 = PDPtr(), tcPDPtr par3 = PDPtr()) const = 0;

  /**
   * Return true, if any of the possible input particle pointer is an
   * exotic quark, e.g. Susy quark; false otherwise.
   */
  virtual bool isExotic(tcPDPtr par1, tcPDPtr par2 = PDPtr(), tcPDPtr par3 = PDPtr()) const = 0;

  //@}

  
  /**
   *  Access the parton weights
   */
   double pwt(long pid) const {
    map<long,double>::const_iterator it = _pwt.find(abs(pid));
    assert( it != _pwt.end() );
    return it->second;
  }

   /**
   * Return true if the two or three particles in input can be the components 
   * of a hadron; false otherwise.
   */
  inline bool canBeHadron(tcPDPtr par1, tcPDPtr par2 , tcPDPtr par3 = PDPtr()) const {
    return (!par3 && canBeMeson(par1,par2)) || canBeBaryon(par1,par2,par3);
  }

  /**
   * Check if can't make a diquark from the partons
   */
  bool canMakeDiQuark(tcPPtr p1, tcPPtr p2) const {
    long id1 = p1->id(), id2 = p2->id();
    return QuarkMatcher::Check(id1) && QuarkMatcher::Check(id2) && id1*id2>0;
  }

  /**
   * Return the particle data of the diquark (anti-diquark) made by the two 
   * quarks (antiquarks) par1, par2.
   * @param par1 (anti-)quark data pointer
   * @param par2 (anti-)quark data pointer
   */
  PDPtr makeDiquark(tcPDPtr par1, tcPDPtr par2) const;

  /**
   * Method to return a pair of hadrons given the PDG codes of
   * two or three constituents
   * @param cluMass The mass of the cluster
   * @param par1 The first constituent
   * @param par2 The second constituent
   * @param par3 The third constituent
   */
  virtual pair<tcPDPtr,tcPDPtr> chooseHadronPair(const Energy cluMass, tcPDPtr par1, 
						 tcPDPtr par2) const;

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
  tcPDPair lightestHadronPair(tcPDPtr ptr1, tcPDPtr ptr2) const;

  /**
   *  Returns the mass of the lightest pair of hadrons with the given particles
   * @param ptr1 is the first  constituent
   * @param ptr2 is the second constituent
   */
  Energy massLightestHadronPair(tcPDPtr ptr1, tcPDPtr ptr2) const  { 
    map<pair<long,long>,tcPDPair>::const_iterator lightest =
      lightestHadrons_.find(make_pair(abs(ptr1->id()),abs(ptr2->id())));
    if(lightest!=lightestHadrons_.end())
      return lightest->second.first->mass()+lightest->second.second->mass();
    else
      return ZERO;
  }

  /**
   * Returns the lightest hadron formed by the given particles.
   *
   * Given the id of two (or three) constituents of a cluster, it returns
   * the  lightest hadron with proper flavour numbers.
   * @param ptr1 is the first  constituent
   * @param ptr2 is the second constituent
   */
   tcPDPtr lightestHadron(tcPDPtr ptr1, tcPDPtr ptr2) const;
  
  /**
   * Return the threshold for a cluster to split into a pair of hadrons.
   * This is normally the mass of the lightest hadron Pair, but can be
   * higher for heavy and exotic clusters
   */
  virtual Energy hadronPairThreshold(tcPDPtr par1, tcPDPtr par2) const=0;

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
  vector<pair<tcPDPtr,double> > hadronsBelowThreshold(Energy threshold,
							      tcPDPtr ptr1, tcPDPtr ptr2) const;

  /**
   * Return the nominal mass of the hadron returned by lightestHadron()
   * @param ptr1 is the first  constituent
   * @param ptr2 is the second constituent 
   * @param ptr3 is the third  constituent 
   */
   Energy massLightestHadron(tcPDPtr ptr1, tcPDPtr ptr2) const {
    // find entry in the table
    pair<long,long> ids(abs(ptr1->id()),abs(ptr2->id()));
    HadronTable::const_iterator tit=_table.find(ids);
    // throw exception if flavours wrong
    if(tit==_table.end()||tit->second.empty())
      throw Exception() <<  "HadronSpectrum::massLightestHadron "
			<< "failed for particle" << ptr1->id()  << " " 
			<< ptr2->id() 
			<< Exception::eventerror;
    // return the mass
    return tit->second.begin()->mass;
  }

  /**
   *  Force baryon/meson selection
   */
  virtual std::tuple<bool,bool,bool> selectBaryon(const Energy cluMass, tcPDPtr par1, tcPDPtr par2) const;

  /**
   *  Returns the mass of the lightest pair of baryons.
   * @param ptr1 is the first  constituent
   * @param ptr2 is the second constituent 
   */
  Energy massLightestBaryonPair(tcPDPtr ptr1, tcPDPtr ptr2) const;
  
  /**
   * Return the weight for the given flavour
   */
   virtual double pwtQuark(const long& id) const = 0;

  virtual double specialQuarkWeight(double quarkWeight, long,
				    const Energy, tcPDPtr, tcPDPtr) const {
    return quarkWeight;
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
  virtual long makeDiquarkID(long id1, long id2, long pspin)  const = 0;

protected:
  /**
   * Return true if the two particles in input can be the components of a meson;
   *false otherwise.
   */
  bool canBeMeson(tcPDPtr par1,tcPDPtr par2)  const;

  /**
   * Return true if the two or three particles in input can be the components 
   * of a baryon; false otherwise.
   */
  virtual bool canBeBaryon(tcPDPtr par1, tcPDPtr par2 , tcPDPtr par3 = PDPtr())  const = 0;


  /**
   *  A sub-function of HadronSpectrum::constructHadronTable().
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
  virtual void constructHadronTable() = 0;

  /**
   * The table of hadron data
   */
  HadronTable _table;

  /**
   *  The PDG codes of the constituent particles allowed
   */
  vector<PDPtr> _partons;

  /**
   *  The PDG codes of the hadrons which cannot be produced in the hadronization
   */
  vector<PDPtr> _forbidden;

  /**
   *  Access to the table of hadrons
   */
   const HadronTable & table() const {
    return _table;
  }
  
  /**
   *  Access to the list of partons
   */
   const vector<PDPtr> & partons() const {
    return _partons;
  }

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
  virtual void insertMeson(HadronInfo a, int flav1, int flav2) = 0;

  /**
   *   Insert a spin\f$\frac12\f$ baryon in the table
   */
  virtual void insertOneHalf(HadronInfo a, int flav1, int flav2);

  /**
   *   Insert a spin\f$\frac32\f$ baryon in the table
   */
  virtual void insertThreeHalf(HadronInfo a, int flav1, int flav2);

  /**
   *  Option for the selection of hadrons below the pair threshold
   */
  unsigned int belowThreshold_;

  /**
   *  The weights for the excited meson multiplets
   */
  vector<vector<vector<double> > > _repwt;

  /**
   * Weights for quarks and diquarks.
   */
  map<long,double> _pwt;

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
   *  Caches of lightest pairs for speed
   */
  //@{
  /**
   * Masses of lightest hadron pair
   */
  map<pair<long,long>,tcPDPair> lightestHadrons_;
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HadronSpectrum & operator=(const HadronSpectrum &) = delete;

};

}

#endif /* Herwig_HadronSpectrum_H */
