// -*- C++ -*-
#ifndef HERWIG_HadronSelector_H
#define HERWIG_HadronSelector_H

#include <ThePEG/Handlers/HandlerBase.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include "CluHadConfig.h"


namespace Herwig {


using namespace ThePEG;


/** \ingroup Hadronization
 *  \class HadronSelector
 *  \brief This class selects the hadron flavours of a cluster decay.
 *  \author Philip Stephens
 *  \author Alberto Ribon
 *
 *  This class performs the selection of a pair of hadrons, 
 *  or in certain cases the selection of a single hadron. These are chosen
 *  with the approriate flavour numbers given by the input.
 * 
 *  This class has implemented three versions of the hadron selection. This
 *  first is given by _ClusterDKMode = 0 and this is the routine implemented
 *  in the Fortran HERWIG program. The second is given by _ClusterDKMode = 1
 *  and this is the method proposed by Kupco. This method should not be
 *  used though as it has been shown to be quite wrong in the prediction of
 *  baryons. The last method, which is the default implementation, is given by
 *  _ClusterDKMode = 2 and this is the new hadronization method developed
 *  for Herwig++.
 */
class HadronSelector: public ThePEG::HandlerBase {

public:

  /**
   * Standard ctors and dtor.
   */
  inline HadronSelector();
  inline HadronSelector(const HadronSelector &);
  virtual ~HadronSelector();

public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

  /**
   * Returns the lightest hadron formed by the given ids.
   *
   * Given the id of two (or three) constituents of a cluster, it returns
   * the id of the lightest hadron with proper flavour numbers.
   * At the moment it does *nothing* in the case that also id3 is present.
   */
  long lightestHadron(const long id1, const long id2, const long id3=0) const;

  /**
   * Return the nominal mass of the hadron with id returned by the previous method.
   */
  Energy massLightestHadron(const long id1, const long id2, const long id3=0) const;

  /**
   * This returns the lightest pair of hadrons given by the flavours.
   *
   * Given the id of two (or three) constituents of a cluster, it returns
   * the id of the two lightest hadrons with proper flavour numbers.
   * Furthermore, the first of the two hadron must have the constituent with
   * id1, and the second must have the constituent with id2. 
   * At the moment it does *nothing* in the case that also id3 is present.
   */
  pair<long,long> 
  lightestHadronPair(const long id1, const long id2, const long id3=0) const;

  Energy massLightestHadronPair(const long id1, const long id2, const long id3=0) 
    const;

  //!< Return the sum of the nominal masses of the two hadrons with id 
  // returned by the previous method.  
  Energy massLightestBaryonPair(const long id1, const long id2) const;

  /** 
   * This method is used to choose a pair of hadrons.
   *
   * Given the mass of a cluster and the ids of its two (or three) 
   * constituents, this returns the pair of ids of the two hadrons with proper
   * flavour numbers. Furthermore, the first of the two hadron must have the 
   * constituent with id1, and the second must have the constituent with id2. 
   * At the moment it does *nothing* in the case that also id3 is present.
   *
   * This routine calls the relevant algorithm depending on what the value of
   * _ClusterDKMode is.
   */
  pair<long,long> chooseHadronPair(const Energy cluMass, const long id1, 
				   const long id2, const long id3=0) 
    throw(Veto, Stop, Exception);

  /**
   * This is the method from HERWIG.
   */
  pair<long,long> hwpp(const Energy, const long, const long, Energy, int);

  /**
   * This is the method proposed by Kupco.
   */
  pair<long,long> kupco(const Energy, const long, const long, Energy, int);

  /**
   * This is the new method developed for Herwig++
   */
  pair<long,long> hw64(const Energy, const long, const long, Energy, int);

  inline double pwtDquark()  const; //!< The down quark weight.
  inline double pwtUquark()  const; //!< The up quark weight.
  inline double pwtSquark()  const; //!< The strange quark weight.
  inline double pwtCquark()  const; //!< The charm quark weight.
  inline double pwtBquark()  const; //!< The bottom quark weight.
  inline double pwtDIquark() const; //!< The diquark weight.

protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;

protected:

  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();

  /**
   * Change all pointers to Interfaced objects to corresponding clones.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return pointers to all Interfaced objects refered to by this.
   */
  inline virtual IVector getReferences();

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<HadronSelector> initHadronSelector;

  /**
   * Private and non-existent assignment operator.
   */
  HadronSelector & operator=(const HadronSelector &);

private:

  /** 
   * Safety check for debugging, usually skipped.
   */
  void safetyCheck(const long id1, const long id2, const long id3=0) const  
    throw(Veto, Stop, Exception);

  /**
   * This method converts the PDG number for the parton into a flavour
   * flavour number used in the maps of this class.
   */
  int  convertIdToFlavour(const long id) const;

  /**
   * Methods converts to the PDG id values for quarks and diquarks
   * from/to the internal consecutive representation of flavours, which starts 
   * from 1 for  d  quark and ends with  bb  diquark at 20.
   * Notice that for quarks d,u,s,c,b the to representation coincides
   * (a part the sign, only positive for the internal flavour representation).
   */
  long convertFlavourToId(const int flavour) const;

  /**
   * This method returns the proper sign ( > 0 hadron; < 0 anti-hadron )
   * for the input PDG id  idHad > 0, suppose to be made by the
   * two constituents of ids: idQ1 and idQ2 (both with proper sign).
   * In the case of failure, it returns 0.
   */
  int  signHadron(const int idQ1, const int idQ2, const int idHad) const;

  /**
   * Build the big table of hadrons and related structures.
   * It calls the method  fillHadronData. 
   */
  void initialize();

  /**
   * This methods is responsible for filling the hadrons input data 
   * (weights). This is the method that one should update when new or updated
   * hadron data is available. Remember always to check whether the
   * new hadrons are present or not in the ThePEG/PDT/EnumParticles.h 
   * if not, then the weights associated to those particles are set
   * automatically to zero, which means that they will never be produced.
   */
  void fillHadronData();
  
private: // data members

  /**
   * Calculates a special weight specific to  a given hadron
   */
  double specialWeight(long); 

  /**
   * Is given id a part of a mixing state
   */
  bool mixingState(long);     

  /**
   * Returns the weight of given mixing state
   */
  double mixingStateWeight(long); 

  /**
   * Parameters
   */
  double _PwtDquark;
  double _PwtUquark;
  double _PwtSquark;
  double _PwtCquark;
  double _PwtBquark;
  double _PwtDIquark;
  double _SngWt; 
  double _DecWt; 

  /** \ingroup Hadronization
   *  \class HadronInfo 
   *  \brief Class used to store all the hadron information for easy access.
   *  \author Philip Stephens
   * 
   *  Notice that: 
   *  - top quark is not present, because its lifetime
   *    is so short to prevent any t-hadrons to be created.
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
   *    an anti-particle is provides in input it is automatically
   *    transform to its particle, simply by taking the module
   *    of its id.
   *  - in order to keep things simpler and to avoid a proliferation 
   *    of classes, although sacrificing a bit of encapsulation,
   *    we have preferred to use below simple classes with all
   *    public members. The main drawback
   *    is that one must be careful to distinguish between the
   *    PDG id (that can be negative and huge in module) with 
   *    the index of vectors. It is useful in this context to use
   *    the methods  convertIdToFlavour  and  covertFlavourToId
   *    where by "flavour" we really mean "index" of related
   *    vectors. 
   */
  class HadronInfo {

    /**
     * Friend method used to print the value of a table element.
     */
    friend PersistentOStream & operator<< ( PersistentOStream & os, 
					    const HadronInfo & hi ) {
      os << hi.id << hi.ptrData << hi.swtef << hi.wt << hi.overallWeight << hi.mass;
      return os;
    }

    /**
     * Friend method used to read in the value of a table element.
     */
    friend PersistentIStream & operator>> ( PersistentIStream & is, HadronInfo & hi ) {
      is >> hi.id >> hi.ptrData >> hi.swtef >> hi.wt >> hi.overallWeight >> hi.mass;
      return is;
    }

  public:

    /**
     * The hadrons id.
     */
    long   id;              

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

    HadronInfo() : id(0), ptrData(tPDPtr()), swtef(1.), wt(1.0), overallWeight(0.0) {}

    bool operator<(const HadronInfo &x) const { return mass < x.mass; }

    void rescale(double x) const { 
      const_cast<HadronInfo*>(this)->overallWeight *= x; 
    }
  };

  /**
   * The type is used to contain all the hadrons info of a given flavour.
   */
  typedef set<HadronInfo> KupcoData;

  /** \ingroup Hadronization
   *  \class Kupco
   *  \brief Class designed to make STL routines easy to use.
   *  \author Philip Stephens
   *
   *  This class is used during the Kupco and Herwig++ routines so that
   *  a list can be generated which is ordered for easy traversal and 
   *  quick access.
   */
  class Kupco {

  public:

    /**
     * id of the quark drawn from the vacuum.
     */
    long idQ;

    /**
     * id of one hadron produced.
     */
    long idHad1;

    /**
     * id of the other hadron produced.
     */
    long idHad2;

    /**
     * Weight factor of this componation.
     */
    double weight;

    /**
     * Simple method to indicate what the < operator means on this class.
     */
    bool operator<(const Kupco &x) const { return weight < x.weight; }

    /**
     * Simple method to indicate what the == operator means on this class.
     */
    bool operator==(const Kupco &x) const { return weight == x.weight; }

    /**
     * Simple method to indicate what the > operator means on this class.
     */
    bool operator>(const Kupco &x) const { return weight > x.weight; }

  };

  /**
   * Make them enums so arrays can be statically allocated
   */

  /**
   * Defines values for array sizes. L,J,N max values for excited mesons.
   */
  enum ConstNumber { NumFlavs = 20,  Lmax = 3, Jmax = 4, Nmax = 4};  

  /**
   * This enum defines the flavour values used in the arrays
   */
  enum Flavs { D = 0, U, S, C, B, DD, DU, UU, DS, US, SS, 
	       DC, UC, SC, CC, DB, UB, SB, CB, BB };  

  /**
   * This object defines the table of hadron data used by the Kupco 
   * and Herwig++ methods.
   */
  KupcoData _table[NumFlavs][NumFlavs];
  
  /**
   * Weights for quarks and diquarks; element [0] empty.
   */
  vector<double> _Pwt;         

  /**
   * Weights for excited mesons;
   */
  double _Repwt[Lmax][Jmax][Nmax];  

  /**
   * whether to use hw64, kupco or hwpp
   */
  int _ClusterDKMode;  

  int _trial;
};
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of HadronsSelector.
 */
template <>
struct BaseClassTrait<Herwig::HadronSelector,1> {
  typedef ThePEG::HandlerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::HadronSelector>: public ClassTraitsBase<Herwig::HadronSelector> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/HadronSelector"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwHadronization.so"; }

};

}

#endif // DOXYGEN

#include "HadronSelector.icc"

#endif /* HERWIG_HadronsSelector_H */
