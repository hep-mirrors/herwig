// -*- C++ -*-
#ifndef ThePEG_BaryonRemnants_H
#define ThePEG_BaryonRemnants_H
//
// This is the declaration of the <!id>BaryonRemnants<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// <!id>BaryonRemnants<!!id> inherits from <!class>RemnantHandler<!!class>
// but can only handle particles without sub-structure with the parton
// density given by a <!class>NoPDF<!!class> object. No , and will
// never give any remnants.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:RemnantHandler.html">RemnantHandler.h</a>,
// <a href="http:NoPDF.html">NoPDF.h</a>.
// 

#include <ThePEG/PDF/RemnantHandler.h>
#include <ThePEG/Utilities/VSelector.h>

namespace Herwig {

using namespace ThePEG;

class BaryonRemnants: public RemnantHandler {

public:

  inline BaryonRemnants();
  inline BaryonRemnants(const BaryonRemnants &);
  virtual ~BaryonRemnants();
  // Standard ctors and dtor

public:

  virtual bool canHandle(tcPDPtr particle,
			 const cPDVector & partons) const;
  // Return true if this remnant handler can handle extracting all
  // specified partons. The BaryonRemnants will return false if any
  // partons are given.

  virtual Lorentz5Momentum generate(PartonBinInstance & pb, const double * r,
				    Energy2 scale,
				    const LorentzMomentum & p) const;
  // Generate the momentum of the extracted parton in the particle cms
  // (but with x still the positive light-cone fraction) as given by
  // the last argument. If the particle is space-like the positive and
  // negative light-cone momenta are sqrt(-m2) and -sqrt(-m2)
  // respectively. If the scale is negative, it means that the doScale
  // in the previous call to nDim() was true, otherwise the given
  // scale should be the virtuality of the extracted parton. Generated
  // quantities which are not returned in the momentum may be saved in
  // the PartonBin for later use. In particular, if the nDim() random
  // numbers are not enough to generate with weight one, the resulting
  // weight should be stored with the remnantWeight() method of the
  // parton bin.

  virtual void createRemnants(PartonBinInstance &) const;
public:

  static void Init();
  // Standard Init function used to initialize the interface.

protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods

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

  static NoPIOClassDescription<BaryonRemnants> initBaryonRemnants;

  BaryonRemnants & operator=(const BaryonRemnants &);
  //  Private and non-existent assignment operator.

};

struct BaryonRemInfo: public RemInfoBase {
  int iq;
  int sign;
  vector<int> flav;
  vector<int> valenceFlav;
  VSelector< pair<int,int> > flavsel;
  bool maybeValence;
};

}

namespace ThePEG {

template <>
struct BaseClassTrait<Herwig::BaryonRemnants,1> {
  typedef RemnantHandler NthBase;
};

template <>
struct ClassTraits<Herwig::BaryonRemnants>
  : public ClassTraitsBase<Herwig::BaryonRemnants> {
  static string className() { return "/Herwig++/BaryonRemnants"; }
  static string library() { return "MRST.so"; }
};

}

#include "BaryonRemnants.icc"

#endif /* ThePEG_BaryonRemnants_H */
