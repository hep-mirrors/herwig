// -*- C++ -*-
#ifndef HERWIG_MECorrection_H
#define HERWIG_MECorrection_H
//
// This is the declaration of the <!id>MECorrection<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This abstract class is the base class from which all concrete <BR>
// matrix element corrections inherit from. It provides a switch for <BR>
// turning on/off the matrix element correction: in this way <BR>
// this switch is automatically available for all the matrix element <BR>
// corrections classes, without need to define it for each new one which is added. <BR> 
// Notice that each object of a class derived from MECorrection must have <BR>
// a reference to: <BR>
// the hard process and hard process plus jet matrix elements; OR <BR>
// the decay process and decay process plus jet matrix elements; <BR>
// according if this ME correction refers to, respectively, <BR>
// a hard <I>2-&GT;N</I> process or a decay <I>1-&GT;N</I> process. <BR>
// 
// This class declares an important pure virtual method: <BR>
//  <!id>softMECorrection.<!!id><BR>
// which must be defined by derived class. Another virtual (but non-pure) <BR>
// method is <!id>hardMECorrection.<!!id><BR>, which does nothing <BR>
// by default, but it could be overriden if Sudakov-suppression of <BR>
// <I>-&GT;N+1<I> matrix elements has to be included.
//
// Notice that:
// <UL>
//  <LI> The approach to Matrix Element correction is quite different than <BR>
//       in Fortran Herwig. In fact, because ThePEG allows to have more <BR>
//       processes competing to each other, we allow the basic <I>-&GT;N</I> <BR>
//       process and the <I>-&GT;N+1</I> one to compete to each other. <BR>
//       Therefore, in the case that the <I>-&GT;N+1</I> process has been <BR>
//       generated, we usually (i.e. without Sudakov suppression of M.E.) <BR>
//       apply the standard showering, without any rejection or soft M.E. <BR>
//       correction. If, instead, the <I>-&GT;N</I> process is generated, <BR>
//       then beside the soft M.E. correction to the hardest emission so far <BR>
//       (exactly as in Fortran Herwig), a further check must be made that <BR>
//       such hardest emission is below some "cut", otherwise the all event <BR>
//       must be rejected (because it will be described by the <I>-&GT;N+1<I>
//       process).         
//  <LI> ***LOOKHERE*** We are assuming that the above "cut", which determines <BR>
//       whether the all generated event will be kept or rejected, is described <BR>
//       in the <!id>KinematicalCut<!!id> object accessible from either <BR> 
//       <!id>hardProcessME()<!!id> or <!id>decayProcessME()<!!id>.  
//  <LI> ***LOOKHERE*** We are assuming that in ThePEG there is a mechanism <BR>
//       (not yet implemented indeed) which allows to know that for each <BR>
//       process with M.E. correction, two processes must compete with each <BR>
//       other: the basic <I>-&GT;N</I> process and the <I>-&GT;N+Jet</I> one.   
// </UL>
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:MECorrections.html">MECorrections.h</a>, <BR>
// <a href="http:DrellYanMECorrection.html">DrellYanMECorrection.h</a>,
// 

#include "ThePEG/Handlers/HandlerBase.h"
#include "ShowerConfig.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/PDT/Decayer.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/CLHEPWrap/Lorentz5Vector.h"


namespace Herwig {

using namespace ThePEG;

class MECorrection: public ThePEG::HandlerBase {

public:

  inline MECorrection();
  inline MECorrection(const MECorrection &);
  virtual ~MECorrection();
  // Standard ctors and dtor.
 
  inline bool isMECorrectionON() const;  
  // It returns true/false if the switch for ME correction is on/off.

  inline MEPtr hardProcessME() const;
  inline MEPtr hardProcessPlusJetME() const;
  // Access the pointers to the matrix elements corresponding
  // to the leading order hard <I>2-&GT;N</I> and <I>2-&GT;N+1</I> processes, 
  // respectively, associated with this ME correction in the case
  // it refers to a hard process; otherwise (that means it refers to a decay 
  // process) those pointers are null.

  inline DecayerPtr decayProcessME() const;
  inline DecayerPtr decayProcessPlusJetME() const;
  // Access the pointers to the matrix elements corresponding
  // to the leading order decay <I>1-&GT;N</I> and <I>1-&GT;N+1</I> processes, 
  // respectively, associated with this ME correction in the case
  // it refers to a decay process; otherwise (that means it refers to a hard
  // process) those pointers are null.

  virtual void hardMECorrection() throw(Veto, Stop, Exception);
  // This method by default does nothing, but it could be overriden
  // by derived classes if Sudakov suppression of hard <I>2-&GT;N+1</I> 
  // or decay <I>1-&GT;N+1</I> processes has to be included,
  // in which case the all event could be rejected.

  virtual bool softMECorrection( Lorentz5Momentum pEmissionHardestSoFar )
    throw(Veto, Stop, Exception) = 0;
  // In the case that the "basic" hard <I>2-&GT;N</I> process or 
  // decay <I>1-&GT;N</I> process has been generated, then for each emission
  // which is the hardest (in <I>Pt</I> w.r.t. to the emitting parent) so far
  // the following two checks are made:
  // --- first, if it is above a certain value (which we assume can be
  //     got from KinematicalCut accessed via the basic process pointer)
  //     then the all event is rejected (because it should be described
  //     by the corresponding <I>-&GT;N+1</I> process);
  // --- second, the soft M.E. correction is applied, which consists in 
  //     accepting/vetoing such emission, which corresponds, respectively, 
  //     to returning false/true in this method (be careful: "true" means
  //     that we have to apply a correction to the normal showering, 
  //     which means, in this case, to reject the emission).
  //     Notice that this method needs only one input, the momentum of the
  //     hardest emission so far, because the other <I>N</I> momenta are 
  //     always the N-outgoing momenta of the hard <I>2-&GT;N</I> process
  //     or decay <I>1-&GT;N</I> process, whose values we can always get 
  //     from the object pointed by <!id>hardProcessME()<!!id> and
  //     <!id>decayProcessME()<!!id>.

public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.

  static void Init();
  // Standard Init function used to initialize the interfaces.

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

  static AbstractClassDescription<MECorrection> initMECorrection;
  // Describe an abstract base class with persistent data.

  MECorrection & operator=(const MECorrection &);
  //  Private and non-existent assignment operator.

  int _correctionMode;
  MEPtr _hardProcessME;
  MEPtr _hardProcessPlusJetME;
  DecayerPtr _decayProcessME;
  DecayerPtr _decayProcessPlusJetME;

};

}

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of MECorrection.
template <>
struct BaseClassTrait<Herwig::MECorrection,1> {
  typedef ThePEG::HandlerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::MECorrection>: public ClassTraitsBase<Herwig::MECorrection> {
  static string className() { return "/Herwig++/MECorrection"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "MECorrection.icc"

#endif /* HERWIG_MECorrection_H */
