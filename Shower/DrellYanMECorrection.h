// -*- C++ -*-
#ifndef HERWIG_DrellYanMECorrection_H
#define HERWIG_DrellYanMECorrection_H
//
// This is the declaration of the <!id>DrellYanMECorrection<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is responsible for the matrix element correction <BR>
// for Drell-Yan processes. <BR>
//
// ***LOOKHERE*** Of course there is no jet real matrix element correction. <BR>
//                This class serves now only just a prototype for the
//                interfaces that all matrix element correctins will have. <BR>
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:MECorrection.html">MECorrection.h</a>.
// 

#include "MECorrection.h"
#include "ShowerConfig.h"
#include "Herwig++/Utilities/GlobalParameters.h"


namespace Herwig {

using namespace ThePEG;

class DrellYanMECorrection: public MECorrection {

public:

  inline DrellYanMECorrection();
  inline DrellYanMECorrection(const DrellYanMECorrection &);
  virtual ~DrellYanMECorrection();
  // Standard ctors and dtor.

  // virtual void hardMECorrection() throw(Veto, Stop, Exception);
  // Override this method only if Sudakov suppression of 
  // hard <I>2-&GT;N+1</I> processes or decay <I>1-&GT;N+1</I> processes 
  // has to be included, in which case the all event could be rejected.

  virtual bool softMECorrection( Lorentz5Momentum pEmissionHardestSoFar )
    throw(Veto, Stop, Exception);
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

  static void Init();
  // Standard Init function used to initialize the interfaces.

protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods.

private:

  static ClassDescription<DrellYanMECorrection> initDrellYanMECorrection;
  // Describe a concrete class with persistent data.

  DrellYanMECorrection & operator=(const DrellYanMECorrection &);
  //  Private and non-existent assignment operator.

};

}

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of DrellYanMECorrection.
template <>
struct BaseClassTrait<Herwig::DrellYanMECorrection,1> {
  typedef Herwig::MECorrection NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::DrellYanMECorrection>: public ClassTraitsBase<Herwig::DrellYanMECorrection> {
  static string className() { return "/Herwig++/DrellYanMECorrection"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "DrellYanMECorrection.icc"

#endif /* HERWIG_DrellYanMECorrection_H */
