// -*- C++ -*-
#ifndef HERWIG_MEee2Z_H
#define HERWIG_MEee2Z_H
//
// This is the declaration of the <!id>MEee2Z<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:.html">.h</a>,
// <a href="http:.html">.h</a>.
// 

#include "ThePEG/MatrixElement/ME2to2Base.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "Herwig++/Models/RSModel/RSModel.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "Herwig++/Helicity/Correlations/ProductionMatrixElement.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
// #include "MEee2Z.fh"
// #include "MEee2Z.xh"

namespace Herwig {
using namespace ThePEG;
using Helicity::ProductionMatrixElement;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using Helicity::VectorWaveFunction;

class MEee2Z: public MEBase {

public:

  inline MEee2Z();
  inline MEee2Z(const MEee2Z &);
  virtual ~MEee2Z();
  // Standard ctors and dtor.

  virtual void constructVertex(tSubProPtr sub);
  // set up the spin correlations

public:

  virtual unsigned int orderInAlphaS() const;
  virtual unsigned int orderInAlphaEW() const;
  // Return the order in respective couplings in which this matrix
  // element is given.

  virtual double me2() const;
  // Return the matrix element for the kinematical configuation
  // previously provided by the last call to setKinematics(), suitably
  // scaled by sHat() to give a dimension-less number.

  virtual Energy2 scale() const;
  // Return the scale associated with the last set phase space point.

  virtual void setKinematics();
  // Set the typed and momenta of the incoming and outgoing partons to
  // be used in subsequent calls to me() and colourGeometries()
  // according to the associated XComb object. If the fun ction is
  // overridden in a sub class the new function must call the base
  // class one first.

  virtual int nDim() const;
  // The number of internal degreed of freedom used in the matrix
  // element. This default version returns 0;

  virtual bool generateKinematics(const double * r);
  // Generate internal degrees of freedom given 'nDim()' uniform
  // random numbers in the interval ]0,1[. To help the phase space
  // generator, the 'dSigHatDR' should be a smooth function of these
  // numbers, although this is not strictly necessary. THe return
  // value should be true of the generation succeeded.

  virtual CrossSection dSigHatDR() const;
  // Return the matrix element squared differential in the variables
  // given by the last call to 'generateKinematics()'.

  virtual void getDiagrams() const;
  // Add all possible diagrams with the add() function.

  inline virtual Selector<DiagramIndex> diagrams(const DiagramVector &) const;
  // With the information previously supplied with the
  // setKinematics(...) method, a derived class may optionally
  // override this method to weight the given diagrams with their
  // (although certainly not physical) relative probabilities.

  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;
  // Return a Selector with possible colour geometries for the selected
  // diagram weighted by their relative probabilities.


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
  
  ProductionMatrixElement HelicityME(Lorentz5Momentum pin[2], Lorentz5Momentum pres,
				     SpinorWaveFunction fin[2],
				     SpinorBarWaveFunction ain[2],
				     VectorWaveFunction vout[3],double&) const;
private:

  static ClassDescription<MEee2Z> initMEee2Z;
  // Describe a concrete class with persistent data.

  MEee2Z & operator=(const MEee2Z &);
  // Private and non-existent assignment operator.

private:

  Ptr<Herwig::Helicity::FFVVertex>::transient_pointer _theFFZVertex;

};

}

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of MEee2Z.
template <>
struct BaseClassTrait<Herwig::MEee2Z,1> {
  typedef MEBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::MEee2Z>
  : public ClassTraitsBase<Herwig::MEee2Z> {
  static string className() { return "Herwig++::MEee2Z"; }
  // Return the class name.
  static string library() { return "HwME.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "MEee2Z.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MEee2Z.tcc"
#endif

#endif /* HERWIG_MEee2Z_H */
