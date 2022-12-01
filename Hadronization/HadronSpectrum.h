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
  
  /** \ingroup Hadronization
   *  \class Kupco
   *  \brief Class designed to make STL routines easy to use.
   *  \author Philip Stephens
   *
   *  This class is used to generate a list of the hadron pairs which can 
   *  be produced that allows easy traversal and quick access.
   */
  class Kupco {
    
  public:
    
    /**
     *  Constructor
     * @param inidQ PDG code of the quark drawn from the vacuum.
     * @param inhad1 ParticleData for the first hadron produced.
     * @param inhad2 ParticleData for the second hadron produced.
     * @param inwgt  The weight for the hadron pair 
     */
     Kupco(tcPDPtr inidQ,tcPDPtr inhad1,tcPDPtr inhad2, Energy inwgt)
      : idQ(inidQ),hadron1(inhad1),hadron2(inhad2),weight(inwgt)
    {}
  
    /**
     * id of the quark drawn from the vacuum.
     */
    tcPDPtr idQ;
    
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
   * Return true if the two or three particles in input can be the components 
   * of a hadron; false otherwise.
   */
  virtual bool canBeHadron(tcPDPtr par1, tcPDPtr par2 , tcPDPtr par3 = PDPtr()) const = 0;

  /**
   *  Check if can't make a diquark from the partons
   */
  virtual bool canMakeDiQuark(tcPPtr p1, tcPPtr p2) const = 0;

  /**
   * Return the particle data of the diquark (anti-diquark) made by the two 
   * quarks (antiquarks) par1, par2.
   * @param par1 (anti-)quark data pointer
   * @param par2 (anti-)quark data pointer
   */
  virtual PDPtr makeDiquark(tcPDPtr par1, tcPDPtr par2) const = 0;

  /**
   * Return the quark flavour which should be considered to set the
   * minimum mass for a minimal cluster splitting.
   */
  virtual long minimalSplitQuark() const = 0;

  /**
   * Method to return a pair of hadrons given the PDG codes of
   * two or three constituents
   * @param cluMass The mass of the cluster
   * @param par1 The first constituent
   * @param par2 The second constituent
   * @param par3 The third constituent
   */
  virtual pair<tcPDPtr,tcPDPtr> chooseHadronPair(const Energy cluMass, tcPDPtr par1, 
						 tcPDPtr par2,tcPDPtr par3 = PDPtr()) const
    = 0;

  /**
   * Select the single hadron for a cluster decay
   * return null pointer if not a single hadron decay
   * @param par1 1st constituent
   * @param par2 2nd constituent
   * @param mass Mass of the cluster
   */
  virtual tcPDPtr chooseSingleHadron(tcPDPtr par1, tcPDPtr par2, Energy mass) const = 0;

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
  virtual pair<tcPDPtr,tcPDPtr> lightestHadronPair(tcPDPtr ptr1, tcPDPtr ptr2,
						   tcPDPtr ptr3 = PDPtr ()) const = 0;

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
   virtual tcPDPtr lightestHadron(tcPDPtr ptr1, tcPDPtr ptr2,
			  tcPDPtr ptr3 = PDPtr ()) const = 0;

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
  virtual vector<pair<tcPDPtr,double> > hadronsBelowThreshold(Energy threshold,
			tcPDPtr ptr1, tcPDPtr ptr2,
			tcPDPtr ptr3 = PDPtr ()) const = 0;

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
      throw Exception() <<  "HadronSpectrum::massLightestHadron "
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
  virtual Energy massLightestBaryonPair(tcPDPtr ptr1, tcPDPtr ptr2) const = 0;
  
  /**
   *  Return the weights for the different quarks and diquarks
   */
  //@{
  /**
   * Return the weight for the given flavour
   */
   virtual double pwtQuark(const long& id) const = 0;

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HadronSpectrum & operator=(const HadronSpectrum &) = delete;

  /**
   * The table of hadron data
   */
  HadronTable _table;

};

}

#endif /* Herwig_HadronSpectrum_H */
