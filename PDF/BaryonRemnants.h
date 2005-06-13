// -*- C++ -*-
#ifndef ThePEG_BaryonRemnants_H
#define ThePEG_BaryonRemnants_H
//
// This is the declaration of the BaryonRemnants class.

#include <ThePEG/PDF/RemnantHandler.h>
#include <ThePEG/Utilities/VSelector.h>

namespace Herwig {

using namespace ThePEG;

/** \ingroup PDF
 *
 *  BaryonRemnants inherits from RemnantHandler but can only handle 
 *  particles without sub-structure with the parton density given by 
 *  a NoPDF object. No, and will never give any remnants.
 *
 *  @see RemnantHandler
 *  @see NoPDF
 */
class BaryonRemnants: public RemnantHandler {

public:

  /**
   * Standard ctors and dtor.
   */
  inline BaryonRemnants();
  inline BaryonRemnants(const BaryonRemnants &);
  virtual ~BaryonRemnants();

public:

  /**
   * Return true if this remnant handler can handle extracting all
   * specified partons. The BaryonRemnants will return false if any
   * partons are given.
   */
  virtual bool canHandle(tcPDPtr particle,
			 const cPDVector & partons) const;

  /**
   * Generate the momentum of the extracted parton in the particle cms
   * (but with x still the positive light-cone fraction) as given by
   * the last argument. If the particle is space-like the positive and
   * negative light-cone momenta are sqrt(-m2) and -sqrt(-m2)
   * respectively. If the scale is negative, it means that the doScale
   * in the previous call to nDim() was true, otherwise the given
   * scale should be the virtuality of the extracted parton. Generated
   * quantities which are not returned in the momentum may be saved in
   * the PartonBin for later use. In particular, if the nDim() random
   * numbers are not enough to generate with weight one, the resulting
   * weight should be stored with the remnantWeight() method of the
   * parton bin.
   */
  virtual Lorentz5Momentum generate(PartonBinInstance & pb, const double * r,
				    Energy2 scale,
				    const LorentzMomentum & p) const;


  // temporary fix here! Need to decide whether to inherit from 
  // ThePEG::BaryonRemnant?
  virtual Lorentz5Momentum generate(PartonBinInstance & pb, const double * r,
				    Energy2 scale, Energy2 ,
				    const LorentzMomentum & parent) const
  {
    return generate(pb,r,scale,parent);
  }
  

  virtual void createRemnants(PartonBinInstance &) const;

public:

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

protected:

  /**
   * Standard clone methods.
   */
  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;

  /**
   * Standard Interfaced virtual functions.
   */
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

  static NoPIOClassDescription<BaryonRemnants> initBaryonRemnants;

  /**
   * Private and non-existent assignment operator.
   */
  BaryonRemnants & operator=(const BaryonRemnants &);

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
