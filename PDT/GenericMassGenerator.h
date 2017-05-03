// -*- C++ -*-
//
// GenericMassGenerator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_GenericMassGenerator_H
#define HERWIG_GenericMassGenerator_H
//
// This is the declaration of the GenericMassGenerator class.
//
#include "ThePEG/PDT/MassGenerator.h" 
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/ParticleData.h"
#include "GenericMassGenerator.fh"
#include "ThePEG/PDT/WidthGenerator.h"
#include "GenericWidthGenerator.fh"
#include "ThePEG/Repository/CurrentGenerator.h"

namespace Herwig {
using namespace ThePEG;

  /**
   *  Declare ModelGenerator class as must be friend to set the particle
   */
  class ModelGenerator;

/** \ingroup PDT
 *
 *  The <code>GenericMassGenerator</code> class is a simple class for the
 *  generation of particle masses in Herwig. It inherits from the 
 *  <code>MassGenerator</code> class of ThePEG and implements a Breit-Wigner
 *  using the width generator to give the running width. 
 *
 *  In general the width generator will be an instance of the
 *  <code>GenericWidthGenerator</code> class which uses the Herwig decayers
 *  based on the <code>DecayIntegrator</code> class to define the shape of the
 *  running width.
 *
 *  This class is designed so that the weight
 *
 *  \f[\int dm^2 \frac{m\Gamma(m)}{(m^2-M^2)^2+m^2\Gamma^2(m)}\f]
 *
 *  can be included in the production of the particle to take off-shell effects into
 *  account. This is the default form of the weight. 
 *  The numerator and running of the width can 
 *  be changed using the BreitWignerShape interface.
 *
 *  @see MassGenerator
 *  @see DecayIntegrator
 *  @see GenericWidthGenerator
 * 
 */
class GenericMassGenerator: public MassGenerator {

  /**
   *  ModelGenerator class as must be friend to set the particle
   */
  friend class ModelGenerator;

public:

  /**
   * Default constructor
   */
  GenericMassGenerator();

  /**
   *  Destructor
   */
  virtual ~GenericMassGenerator();

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

public:

  /**
   * Return true if this mass generator can handle the given particle type.
   * @param part The particle data pointer of the particle.
   * @return True ig this class can handle the particle and false otherwise
   */
  bool accept(const ParticleData & part) const;

  /** @name Members to generate the mass of a particle instance */
  //@{
  /**
   * Generate a mass using the default limits.
   * @param part The particle data pointer of the particle.
   * @return The mass of the particle instance.
   */
  Energy mass(const ParticleData & part) const {
    return mass(part,lowerMass_,upperMass_);
  }

  /**
   * Generate a mass using specified limits
   * @param part The particle data pointer of the particle.
   * @param low The lower limit on the particle's mass.
   * @param upp The upper limit on the particle's mass.
   * @return The mass of the particle instance.
   */
  Energy mass(const ParticleData & part,
	      const Energy low,const Energy upp) const {
    if(upp<low) return low;
    Energy output;
    int ntry=0; double wgt=0.;
    do {
      ++ntry;
      output=mass(wgt,part,low,upp,3);
      if(wgt>maxWgt_) maxWgt_=wgt;
    }
    while(maxWgt_*(UseRandom::rnd())>wgt&&ntry<nGenerate_);
    return (ntry>=nGenerate_) ? mass_ : output;
  }

  /**
   * Return a mass with the weight using the default limits.
   * @param part The particle data pointer of the particle.
   * @param wgt The weight for this mass.
   * @param r   The random number used for the weight
   * @return The mass of the particle instance.
   */
  Energy mass(double & wgt, const ParticleData & part, 
	      double r=UseRandom::rnd()) const {
    return mass(wgt,part,lowerMass_,upperMass_,r);
  }
  
  /**
   * Return a mass with the weight using the specified limits.
   * @param part The particle data pointer of the particle.
   * @param low The lower limit on the particle's mass.
   * @param upp The upper limit on the particle's mass.
   * @param wgt The weight for this mass.
   * @param r   The random number used for the weight
   * @return The mass of the particle instance.
   */
  Energy mass(double & wgt, const ParticleData & part,
	      const Energy low,const Energy upp,
	      double r=UseRandom::rnd()) const {
    return mass(wgt,part,low,upp,BWShape_,r);
  }
  
  /**
   * Weight for the factor.
   * @param q The mass of the instance
   * @return The weight.
   */
  virtual double weight(Energy q) const {
    return weight(q,BWShape_);
  }

  /**
   *  Return the full weight
   */
  virtual InvEnergy2 BreitWignerWeight(Energy q) {
    return BreitWignerWeight(q,BWShape_);
  }
  //@}

  /**
   * Output the initialisation info for the database
   */
  virtual void dataBaseOutput(ofstream &,bool);

public:

  /** @name Access to particle properties */
  //@{
  /**
   * The running width.
   * @param q The mass for the calculation of the running width
   * @return The running width.
   */
  pair<Energy,Energy> width(Energy q,int shape) const;

  /**
   * Lower limit on the mass
   */
  Energy lowerLimit() const {return lowerMass_;}

  /**
   * Upper limit on the mass
   */
  Energy upperLimit() const {return upperMass_;}

  /**
   * Default mass
   */
  Energy nominalMass() const {return mass_;}

  /**
   * Default Width
   */
  Energy nominalWidth() const {return width_;}

protected:

  /**
   * Return a mass with the weight using the specified limits.
   * @param low The lower limit on the particle's mass.
   * @param upp The upper limit on the particle's mass.
   * @param wgt The weight for this mass.
   * @param shape The type of shape to use
   * @param r   The random number used for the weight
   * @return The mass of the particle instance.
   */
  virtual Energy mass(double & wgt, const ParticleData & ,
		      const Energy low,const Energy upp, int shape,
		      double r=UseRandom::rnd()) const {
    // calculate the mass square using fixed width BW
    Energy  lo=max(low,lowerMass_),up=min(upp,upperMass_);
    double  rhomin=atan2((lo*lo-mass2_),mWidth_);
    double  rhomax=atan2((up*up-mass2_),mWidth_)-rhomin;
    double  rho=rhomin+rhomax*r;
    Energy2 q2 = mass2_+mWidth_*tan(rho);
    Energy  q = sqrt(q2);  
    wgt = rhomax*weight(q,shape);
    // return the mass
    return q;
  }

  /**
   * Return a mass with the weight using the default limits.
   * @param part The particle data pointer of the particle.
   * @param wgt The weight for this mass.
   * @param shape The type of shape to use
   * @param r   The random number used for the weight
   * @return The mass of the particle instance.
   */
  Energy mass(double & wgt, const ParticleData & part, int shape,
	      double r=UseRandom::rnd()) const {
    return mass(wgt,part,lowerMass_,upperMass_,shape,r);
  }

  /**
   * Weight for the factor.
   * @param q The mass of the instance
   * @param shape The type of shape to use as for the BreitWignerShape interface
   * @return The weight.
   */
  inline virtual double weight(Energy q, int shape) const {
    Energy2 q2 = q*q;
    Energy4 sq=sqr(q2-mass2_);
    pair<Energy,Energy> gam=width(q,shape);
    // finish the calculation of the width
    Energy2 num;
    if(shape==2)      num = mass_*gam.first;
    else if(shape==3) num = mass_*gam.first;
    else              num = q    *gam.first;
    Energy4 den = (shape==2) ? sq+mass2_*gam.second*gam.second : sq+q2*gam.second*gam.second;
    return num/den*(sq+mWidth_*mWidth_)/Constants::pi/mWidth_;
  }

  /**
   *  Return the full weight
   */
  virtual InvEnergy2 BreitWignerWeight(Energy q, int shape) const {
    Energy2 q2 = q*q;
    Energy4 sq=sqr(q2-mass2_);
    pair<Energy,Energy> gam=width(q,shape);
    // finish the calculation of the width
    Energy2 num;
    if(shape==2)      num = mass_*gam.first;
    else if(shape==3) num = mass_*gam.first;
    else              num = q    *gam.first;
    Energy4 den = (shape==2) ? sq+mass2_*gam.second*gam.second : sq+q2*gam.second*gam.second;
    return num/den/Constants::pi;
  }

  /**
   *  Accesss to the particle
   */
  tcPDPtr particle() const {return particle_;}

  /**
   * Set the particle
   */
  void particle(tPDPtr in) {particle_ = in;}

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans)
   ;

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<GenericMassGenerator> initGenericMassGenerator;

  /**
   * Private and non-existent assignment operator.
   */
  GenericMassGenerator & operator=(const GenericMassGenerator &);

private:
 
  /**
   * Helper function for the interface
   */
  void setParticle(string);

  /**
   * Helper function for the interface
   */
  string getParticle() const;

private:

  /**
   * The maximum weight for unweighting when generating the mass.
   */
  mutable double maxWgt_;

  /**
   * parameter controlling the shape of the Breit-Wigner
   */
  int BWShape_;

  /**
   * Number of attempts to generate the mass.
   */
  int nGenerate_;

private:

  /**
   * Pointer to the particle
   */
  tPDPtr particle_;

  /**
   * Lower limit on the particle's mass
   */
  Energy lowerMass_;

  /**
   * Upper limit on the particle's mass
   */
  Energy upperMass_;

  /**
   * Mass of the particle
   */
  Energy mass_;

  /**
   * Width of the particle
   */
  Energy width_; 

  /**
   * Mass of the particle squared.
   */
  Energy2 mass2_;

  /**
   * Mass of the particle times the width.
   */
  Energy2 mWidth_;

  /**
   * Number of weights to generate when initializing
   */
  int nInitial_;

  /**
   * Whether or not to initialize the GenericMassGenerator
   */
  bool initialize_;

  /**
   * Whether or not to output the data to a file
   */
  bool output_;

  /**
   * Pointer to the width generator
   */
  WidthGeneratorPtr widthGen_;

  /**
   * Pointer to the width generator
   */
  GenericWidthGeneratorPtr widthGenB_;

  /**
   *  Option for the treatment of the width
   */
  bool widthOpt_;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of GenericMassGenerator.
 */
template <>
struct BaseClassTrait<Herwig::GenericMassGenerator,1> {
  /** Typedef of the base class of GenericMassGenerator. */
  typedef MassGenerator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::GenericMassGenerator>
  : public ClassTraitsBase<Herwig::GenericMassGenerator> {
  /** Return the class name. */
  static string className() { return "Herwig::GenericMassGenerator"; }
};

/** @endcond */

}

#endif /* HERWIG_GenericMassGenerator_H */
