// -*- C++ -*-
#ifndef ThePEG_MEqq2gZ2ll_H
#define ThePEG_MEqq2gZ2ll_H
//
// This is the declaration of the <!id>MEqq2gZ2ll<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// The <!id>MEqq2gZ2ll<!!id> class implements the
// <i>e<sup>+</sup>e<sup>-</sup> -&gt; gamma/Z<sup>0</sup> -&gt;
// q+qbar</i> matrix element. Both the continuum and
// <i>Z<sup>0</sup></i> pole term as well as the interference term is
// included. Although not a strict QCD matrix element the cass
// inherits from <!class>ME2to2Base<!!class>, mainly to inherit the
// parameter for the number of active quark flavours.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:ME2to2QCD.html">ME2to2QCD.h</a>.
// 

#include "ThePEG/MatrixElement/ME2to2QCD.h"

using namespace ThePEG;

namespace Herwig {

class MEqq2gZ2ll: public ME2to2QCD {

public:

  MEqq2gZ2ll();
  MEqq2gZ2ll(const MEqq2gZ2ll &);
  virtual ~MEqq2gZ2ll();
  // Standard ctors and dtor

public:

  virtual unsigned int orderInAlphaS() const;
  virtual unsigned int orderInAlphaEW() const;
  // Return the order in respective couplings in which this matrix
  // element is given. Returns 0 and 2 respectively.

  virtual double me2() const;
  // Return the matrix element for the kinematical configuation
  // previously provided by the last call to setKinematics().

  virtual void getDiagrams() const;
  // Add all possible diagrams with the add() function.

  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;
  // Return a Selector with possible colour geometries for the selected
  // diagram weighted by their relative probabilities.

  virtual Selector<DiagramIndex> diagrams(const DiagramVector &) const;
  // Weight the given diagrams with their relative probabilities.

  virtual Energy2 scale() const;
  // Return the scale associated with the last set phase space point.

public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.

  static void Init();
  // Standard Init function used to initialize the interface.

protected:

  virtual IBPtr clone() const;
  virtual IBPtr fullclone() const;
  // Standard clone methods

  virtual void doinit() throw(InitException);
  // Standard Interfaced virtual functions.

protected:

  vector<double> coefs;
  // Constants for the different terms set from the StandardModel in
  // the init() function.

  Energy2 mZ2;
  Energy2 GZ2;
  // The mass squared and width squared of the Z0.

  mutable double lastCont;
  mutable double lastBW;
  // The last continuum and Breit-Wigner terms to be used to select
  // primary diagram.

private:

  static ClassDescription<MEqq2gZ2ll> initMEqq2gZ2ll;

  MEqq2gZ2ll & operator=(const MEqq2gZ2ll &);
  //  Private and non-existent assignment operator.

};

}

namespace ThePEG {

template <>
struct BaseClassTrait<Herwig::MEqq2gZ2ll,1> {
  typedef ME2to2QCD NthBase;
};

template <>
struct ClassTraits<Herwig::MEqq2gZ2ll>
  : public ClassTraitsBase<Herwig::MEqq2gZ2ll> {
  static string className() { return "/Herwig++/MEqq2gZ2ll"; }
  static string library() { return "libHwME.so"; }
};

}

#include "MEqq2gZ2ll.icc"

#endif /* ThePEG_MEqq2gZ2ll_H */
