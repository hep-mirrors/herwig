// -*- C++ -*-
#ifndef HERWIG_GenericMassGenerator_H
#define HERWIG_GenericMassGenerator_H
//
// This is the declaration of the GenericMassGenerator class.
//
#include "ThePEG/PDT/MassGenerator.h" 
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/WidthGenerator.h"
#include "GenericMassGenerator.fh"
#include "ThePEG/Repository/CurrentGenerator.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup PDT
 *
 *  The <code>GenericMassGenerator</code> class is a simple class for the
 *  generation of particle masses in Herwig++. It inherits from the 
 *  <code>MassGenerator</code> class of ThePEG and implements a Breit-Wigner
 *  using the width generator to give the running width. 
 *
 *  In general the width generator will be an instance of the
 *  <code>GenericWidthGenerator</code> class which uses the Herwig++ decayers
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

public:

  /**
   * Default constructor
   */
  inline GenericMassGenerator();

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
  inline Energy mass(const ParticleData & part) const;

  /**
   * Generate a mass using specified limits
   * @param part The particle data pointer of the particle.
   * @param low The lower limit on the particle's mass.
   * @param upp The upper limit on the particle's mass.
   * @return The mass of the particle instance.
   */
  inline Energy mass(const ParticleData & part,const Energy low,const Energy upp) const;

  /**
   * Return a mass with the weight using the default limits.
   * @param part The particle data pointer of the particle.
   * @param wgt The weight for this mass.
   * @return The mass of the particle instance.
   */
  inline Energy mass(const ParticleData & part,double & wgt) const;

  /**
   * Return a mass with the weight using the specified limits.
   * @param part The particle data pointer of the particle.
   * @param low The lower limit on the particle's mass.
   * @param upp The upper limit on the particle's mass.
   * @param wgt The weight for this mass.
   * @return The mass of the particle instance.
   */
  inline Energy mass(const ParticleData & part,double & wgt, 
		     const Energy low,const Energy upp) const;

  /**
   * Weight for the factor.
   * @param mass The mass of the instance
   * @return The weight.
   */
  inline virtual double weight(Energy mass) const;
  //@}

  /**
   * Output the initialisation info for the database
   */
  void dataBaseOutput(ofstream &);

public:

  /** @name Access to particle properties */
  //@{
  /**
   * The running width.
   * @param mass The mass for the calculation of the running width
   * @return The running width.
   */
  inline Energy width(Energy mass) const;

  /**
   * Lower limit on the mass
   */
  inline Energy lowerLimit() const;

  /**
   * Upper limit on the mass
   */
  inline Energy upperLimit() const;

  /**
   * Default mass
   */
  inline Energy nominalMass() const;

  /**
   * Default Width
   */
  inline Energy nominalWidth() const;

protected:

  /**
   * Return a mass with the weight using the specified limits.
   * @param part The particle data pointer of the particle.
   * @param low The lower limit on the particle's mass.
   * @param upp The upper limit on the particle's mass.
   * @param wgt The weight for this mass.
   * @param shape The type of shape to use
   * @return The mass of the particle instance.
   */
  inline Energy mass(const ParticleData & part,double & wgt, 
		     const Energy low,const Energy upp, int shape) const;

  /**
   * Return a mass with the weight using the default limits.
   * @param part The particle data pointer of the particle.
   * @param wgt The weight for this mass.
   * @param shape The type of shape to use
   * @return The mass of the particle instance.
   */
  inline Energy mass(const ParticleData & part,double & wgt,int shape) const;

  /**
   * Weight for the factor.
   * @param mass The mass of the instance
   * @param shape The type of shape to use as for the BreitWignerShape interface
   * @return The weight.
   */
  inline virtual double weight(Energy mass,int shape) const;

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
  virtual void doinit() throw(InitException);
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

 protected:

  /**
   * The maximum weight for unweighting when generating the mass.
   */
  mutable double _maxwgt;

  /**
   * parameter controlling the shape of the Breit-Wigner
   */
  int _BWshape;

  /**
   * Number of attempts to generate the mass.
   */
  int _ngenerate;

private:
 
  /**
   * Pointer to the particle
   */
  PDPtr _particle;

  /**
   * Pointer to the anti-particle
   */
  PDPtr _antiparticle;

  /**
   * Lower limit on the particle's mass
   */
  Energy _lowermass;

  /**
   * Upper limit on the particle's mass
   */
  Energy _uppermass;

  /**
   * Mass of the particle
   */
  Energy _mass;

  /**
   * Width of the particle
   */
  Energy _width; 

  /**
   * Mass of the particle squared.
   */
  Energy2 _mass2;

  /**
   * Mass of the particle times the width.
   */
  Energy2 _mwidth;

  /**
   * Number of weights to generate when initializing
   */
  int _ninitial;

  /**
   * Whether or not to initialize the GenericMassGenerator
   */
  bool _initialize;

  /**
   * Pointer to the width generator
   */
  WidthGeneratorPtr _widthgen;
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

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
  static string className() { return "Herwig++::GenericMassGenerator"; }
};

}

#include "GenericMassGenerator.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "GenericMassGenerator.tcc"
#endif

#endif /* HERWIG_GenericMassGenerator_H */
