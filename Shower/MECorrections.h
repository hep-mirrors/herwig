// -*- C++ -*-
#ifndef HERWIG_MECorrections_H
#define HERWIG_MECorrections_H
//
// This is the declaration of the <!id>MECorrections<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is responsible for managing the collection of <BR>
// matrix element corrections. It also has two switches for <BR>
// turning on/off all the corrections, and for allow/disallow <BR>
// the compositions of such corrections (like for example <BR>
// the top production correction and the top decay correction <BR>
// in the same event).
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:MECorrection.html">MECorrection.h</a>.
// 

#include "Pythia7/Handlers/HandlerBase.h"
#include "ShowerConfig.h"
#include "Herwig++/Config/GlobalParameters.h"
#include "MECorrection.h"
#include "Pythia7/MatrixElement/MEBase.h"
#include "Pythia7/PDT/Decayer.h"


namespace Herwig {

using namespace Pythia7;

class MECorrections: public Pythia7::HandlerBase {

public:

  inline MECorrections();
  inline MECorrections(const MECorrections &);
  virtual ~MECorrections();
  // Standard ctors and dtor.

  inline bool isMECorrectionsON() const;  
  // It returns true/false if the overall switch for ME corrections is on/off.

  inline bool isComposeMECorrectionsON() const;  
  // It returns true if <!id>isMECorrectionsON()<!!id> is true and the composition 
  // for ME corrections switch is on; false otherwise.
  // For example, suppose that both the ME corrections for top production and
  // top decay are present and with switch on; if <!id>isComposeMECorrectionsON()<!!id>
  // is true, then in <I>t tbar</I> events we have to apply the ME correction at the
  // production and also (and independently) the ME correction for top decay.
  // Of course, independence between the two (or more) kinds of ME corrections
  // is always assumed.

  tMECorrectionPtr getMECorrection( const tMEPtr hardProcess ) const;
  tMECorrectionPtr getMECorrection( const tDecayerPtr decayProcess ) const;
  // Given the id of a certain process, hard subprocess or decay respectively, 
  // it returns the pointer to the matrix element correction for such process, 
  // if defined and switched on; otherwise a null pointer is returned.

public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.

  static void Init();
  // Standard Init function used to initialize the interfaces.

protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods.

protected:

  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void doinitrun();
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.

  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.

private:

  static ClassDescription<MECorrections> initMECorrections;
  // Describe a concrete class with persistent data.

  MECorrections & operator=(const MECorrections &);
  //  Private and non-existent assignment operator.

  void initializeRun();
  // It is called once, at the beginning of the run, to build the map of 
  // processes (hard subprocesses and decays) that have matrix element corrections. 

  vector<MECorrectionPtr> _vecMECorrection;  // collection of (pointers to) ME corrections
  int _MECorrectionsMode;
  int _composeMECorrectionsMode;

  map<tMEPtr,int> _mapHardProcesses; 
  map<tMEPtr,int> _mapHardPlusJetProcesses; 
  // key = hard process ME; value = index position in _vecMECorrection
  // Notice that two maps are necessary: one for the basic hard <I>2-&GT;N</I>
  // process, and one for the <I>2-&GT;N+1</I> process, because these two
  // processes compete with each other, and events can be produced with one
  // or the other. Notice also that the same process could appear in both
  // maps (for example in the case of simultaneous e+ e- -&GT; 2,3,4,... jets
  // matching with parton showers, in Catani-Krauss-Webber paper).

  map<tDecayerPtr,int> _mapDecayProcesses; 
  map<tDecayerPtr,int> _mapDecayPlusJetProcesses; 
  // key = process ME; value = index position in _vecMECorrection
  // A similar comment, as above, applies here for the need of two maps.

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of MECorrections.
template <>
struct BaseClassTrait<Herwig::MECorrections,1> {
  typedef Pythia7::HandlerBase NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::MECorrections>: public ClassTraitsBase<Herwig::MECorrections> {
  static string className() { return "/Herwig++/MECorrections"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "MECorrections.icc"

#endif /* HERWIG_MECorrections_H */
