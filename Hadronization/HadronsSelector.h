// -*- C++ -*-
#ifndef HERWIG_HadronsSelector_H
#define HERWIG_HadronsSelector_H
//
// This is the declaration of the <!id>HadronsSelector<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class performs the selection of, normally, a pair of hadrons, <BR>
// or, exceptionally, of a single hadron, with the right flavour numbers. <BR>
// The choice of a pair of hadrons is made using Kupco's method. <BR>
//

#include "Pythia7/Handlers/HandlerBase.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"
#include "CluHadConfig.h"


namespace Herwig {


using namespace Pythia7;


class HadronsSelector: public Pythia7::HandlerBase {

public:

  inline HadronsSelector();
  inline HadronsSelector(const HadronsSelector &);
  virtual ~HadronsSelector();
  // Standard ctors and dtor.

public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.

  static void Init();
  // Standard Init function used to initialize the interfaces.

  long lightestHadron(const long id1, const long id2, const long id3=0) const;
  // Given the id of two (or three) constituents of a cluster, it returns
  // the id of the lightest hadron with proper flavour numbers.
  // At the moment it does *nothing* in the case that also id3 is present.

  Energy massLightestHadron(const long id1, const long id2, const long id3=0) const;
  // Return the nominal mass of the hadron with id returned by the previous method.

  pair<long,long> 
  lightestHadronsPair(const long id1, const long id2, const long id3=0) const;
  // Given the id of two (or three) constituents of a cluster, it returns
  // the id of the two lightest hadrons with proper flavour numbers.
  // Furthermore, the first of the two hadron must have the constituent with
  // id1, and the second must have the constituent with id2. 
  // At the moment it does *nothing* in the case that also id3 is present.

  Energy massLightestHadronsPair(const long id1, const long id2, const long id3=0) 
    const;
  // Return the sum of the nominal masses of the two hadrons with id returned 
  // by the previous method.

  pair<long,long> chooseHadronsPair(const Energy cluMass, const long id1, const long id2, 
				    const long id3=0) throw(Veto, Stop, Exception);
  // Given the mass of a cluster and the ids of its two (or three) constituents, 
  // it returns the pair of ids of the two hadrons with proper flavour numbers.
  // Furthermore, the first of the two hadron must have the constituent with
  // id1, and the second must have the constituent with id2. 
  // At the moment it does *nothing* in the case that also id3 is present.

  inline double pwtDquark()  const;
  inline double pwtUquark()  const;
  inline double pwtSquark()  const;
  inline double pwtCquark()  const;
  inline double pwtBquark()  const;
  inline double pwtDIquark() const;
  // Simple "get" methods to provide the weight parameters to other classes.
  // Probably only the class <!class>ClusterFissioner<!!class> uses the above 
  // methods (indeed it uses the first three only)

protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods.

protected:

  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.

  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.

private:

  static ClassDescription<HadronsSelector> initHadronsSelector;
  // Describe a concrete class with persistent data.

  HadronsSelector & operator=(const HadronsSelector &);
  //  Private and non-existent assignment operator.

private:

  void safetyCheck(const long id1, const long id2, const long id3=0) const  
    throw(Veto, Stop, Exception);
  // Safety check for debugging, usually skipped.

  int  convertIdToFlavour(const long id) const;
  long convertFlavourToId(const int flavour) const;
  // Methods to convert to/from the PDG id values for quarks and diquarks
  // from/to the internal consecutive representation of flavours, which starts 
  // from 1 for  d  quark and ends with  bb  diquark at 20.
  // Notice that for quarks d,u,s,c,b the to representation coincides
  // (a part the sign, only positive for the internal flavour representation).

  int  signHadron(const int idQ1, const int idQ2, const int idHad) const;
  // This method returns the proper sign ( > 0 hadron; < 0 anti-hadron )
  // for the input PDG id  idHad > 0, suppose to be made by the
  // two constituents of ids: idQ1 and idQ2 (both with proper sign).
  // In the case of failure, it returns 0.

  void initialize();
  // Build the big table of hadrons and related structures.
  // It calls the method  fillDataHadrons. 

  int fillDataHadrons();
  // Methods responsibles to fill the hadrons input data (weights).
  // This is the method that one should update when new or updated
  // hadron data is available. Remember always to check whether the
  // new hadrons are present or not in the Pythia7/PDT/EnumParticles.h 
  // if not, then the weights associated to those particles are set
  // automatically to zero, which means that they will never be produced.
  
private: // data members

  // Parameters
  double _PwtDquark;
  double _PwtUquark;
  double _PwtSquark;
  double _PwtCquark;
  double _PwtBquark;
  double _PwtDIquark;
  double _SngWt; 
  double _DecWt; 

  // Types and objects used to implement Kupco's method for
  // selecting the pair of hadrons with proper flavours.
  // Notice that: 
  // --- top quark is not present, because its lifetime
  //     is so short to prevent any t-hadrons to be created.
  // --- the big vector of hadron data, _vecHad, has the first
  //     [0] element empty; similarly, the location table
  //     _locHad has its first element [0][0] empty, as well
  //     many of its elements (the ones that would correspond
  //     to unphysical  diquark-diquark  bound state).
  //     Differently from _vecHad and _locHad, Kupco's table 
  //     kupcoTable is instead filled starting from the first 
  //     element [0] (simply because the index manipulation 
  //     arguments for _vecHad and _locHad do not apply for 
  //     kupcoTable).
  // --- the hadrons in _vecHad can be filled in any ordered 
  //     w.r.t. the mass value, and flavours for different
  //     groups (for instance, (u,s) hadrons don't need to
  //     be placed after (d,s) or any other flavour), but 
  //     all hadrons with the same flavours must be consecutive
  //     ( for instance you cannot alternate hadrons of type
  //       (d,s) with those of flavour (u,s) ).
  //     Furthermore, it is assumed that particle and antiparticle
  //     have the same weights, and therefore only one of them
  //     must be entered in _vecHad: we have chosen to refer
  //     to the particle, defined as PDG id > 0, although if
  //     an anti-particle is provides in input it is automatically
  //     transform to its particle, simply by taking the module
  //     of its id.
  // --- in order to keep things simpler and to avoid a proliferation 
  //     of classes, although sacrificing a bit of encapsulation,
  //     we have preferred to use below simple classes with all
  //     public members. The main drawback
  //     is that one must be careful to distinguish between the
  //     PDG id (that can be negative and huge in module) with 
  //     the index of vectors. It is useful in this context to use
  //     the methods  convertIdToFlavour  and  covertFlavourToId
  //     where by "flavour" we really mean "index" of related
  //     vectors. 

  class HadronInfo {
    friend PersistentOStream & operator<< ( PersistentOStream & os, const HadronInfo & hi ) {
      os << hi.id << hi.ptrData << hi.swtef << hi.wt << hi.overallWeight << hi.mass;
      return os;
    }
    friend PersistentIStream & operator>> ( PersistentIStream & is, HadronInfo & hi ) {
      is >> hi.id >> hi.ptrData >> hi.swtef >> hi.wt >> hi.overallWeight >> hi.mass;
      return is;
    }
  public:
    long   id;              
    tPDPtr ptrData;        // pointer to ParticleData, to get the spin, etc...
    double swtef;          // singlet/decuplet/orbital factor 
    double wt;             // mixing factor
    double overallWeight;  // (2*J+1)*wt*swtef
    Energy mass;
  };
  enum ConstNumber { NumHadrons = 4000, NumQuarksDiquarks = 20, MaxNumTries = 1000,
                     Lmax = 3, Jmax = 4, Nmax = 4};  // L,J,N max values for excited mesons.
  vector<HadronInfo> _vecHad;  // Only hadrons with id > 0; no anti-hadrons!
                               // The first element _vecHad[0] is empty.
  vector<double> _Pwt;         // Weights for quarks and diquarks; element [0] empty.
  enum Flavours { D=1, U, S, C, B, DD, UD, UU, SD, SU, SS,
                  CD, CU, CS, CC, BD, BU, BS, BC, BB };
  class HadronLocation {
    friend PersistentOStream & operator<< ( PersistentOStream & os, const HadronLocation & hl ) {
      os << hl.first << hl.last << hl.lightest;
      return os;
    }
    friend PersistentIStream & operator>> ( PersistentIStream & is, HadronLocation & hl ) {
      is >> hl.first >> hl.last >> hl.lightest;
      return is;
    }
  public:
    int first;     // position in vecRes of the first    resonance of that flavour
    int last;      // position in vecRes of the last     resonance of that flavour
    int lightest;  // position in vecRes of the lightest resonance of that flavour
  };
  vector< vector<HadronLocation> > _locHad;
  // Only the following elements of  _locHad[i][j]  are really used:
  //   D <= i <= B  and   i <= j <= NumQuarksDiquarks
  // and symmetric _locHad[j][i]  

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of HadronsSelector.
template <>
struct BaseClassTrait<Herwig::HadronsSelector,1> {
  typedef Pythia7::HandlerBase NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::HadronsSelector>: public ClassTraitsBase<Herwig::HadronsSelector> {
  static string className() { return "/Herwig++/HadronsSelector"; }
  // Return the class name.
  static string library() { return "libHwHadronization.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "HadronsSelector.icc"

#endif /* HERWIG_HadronsSelector_H */
