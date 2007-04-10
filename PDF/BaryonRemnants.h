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

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline BaryonRemnants();

  /**
   * Destructor.
   */
  virtual ~BaryonRemnants();
  //@}

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


  /**
   * temporary fix here! Need to decide whether to inherit from 
   * ThePEG::BaryonRemnant?
   */
  virtual Lorentz5Momentum generate(PartonBinInstance & pb, const double * r,
				    Energy2 scale, Energy2 ,
				    const LorentzMomentum & parent) const
  {return generate(pb,r,scale,parent);}
  
  /**
   *  Boost the remnants
   */
  virtual void boostRemnants(PartonBinInstance &) const;

public:

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:

  /**
   *  describe a class without persistent data
   */
  static NoPIOClassDescription<BaryonRemnants> initBaryonRemnants;

  /**
   * Private and non-existent assignment operator.
   */
  BaryonRemnants & operator=(const BaryonRemnants &);
};
}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of BaryonRemnants.
 */
template <>
struct BaseClassTrait<Herwig::BaryonRemnants,1> {
  /** Typedef of the base class of BaryonRemnants. */
  typedef RemnantHandler NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::BaryonRemnants>
  : public ClassTraitsBase<Herwig::BaryonRemnants> {
   /** Return the class name. */
  static string className() { return "Herwig++::BaryonRemnants"; }
   /**
    * Return the name of the shared library to be loaded to get
    * access to this class and every other class it uses
    * (except the base class).
    */
  static string library() { return "HwMRST.so"; }
};

/** @endcond */

}

#include "BaryonRemnants.icc"

#endif /* ThePEG_BaryonRemnants_H */
