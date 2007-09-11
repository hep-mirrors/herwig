// -*- C++ -*-
#ifndef HERWIG_SMWZDecayer_H
#define HERWIG_SMWZDecayer_H
//
// This is the declaration of the SMWZDecayer class.
//
#include "Herwig++/Decay/DecayIntegrator.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "SMWZDecayer.fh"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/** \ingroup Decay
 *
 *  The <code>SMWZDecayer</code> is designed to perform the decay of the 
 *  W and Z bosons to the Standard Model fermions. In principle it can also
 *  be used for these decays in any model.
 *
 * @see DecayIntegrator
 * 
 */
class SMWZDecayer: public DecayIntegrator {

public:

  /**
   * Default constructor.
   */
  SMWZDecayer();

  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const PDVector & children) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param vertex Output the information on the vertex for spin correlations
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(bool vertex, const int ichan, const Particle & part,
		     const ParticleVector & decay) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

public:

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

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<SMWZDecayer> initSMWZDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  SMWZDecayer & operator=(const SMWZDecayer &);

 private:

  /**
   * Pointer to the Z vertex
   */
  FFVVertexPtr _zvertex;

  /**
   * Pointer to the W vertex
   */
  FFVVertexPtr _wvertex;

  /**
   * maximum weights for the different integrations
   */
  //@{
  /**
   *  Weights for the Z to quarks decays.
   */
  vector<double> _zquarkwgt;

  /**
   *  Weights for the Z to leptons decays.
   */
  vector<double> _zleptonwgt;

  /**
   *  Weights for the W to quarks decays.
   */
  vector<double> _wquarkwgt;

  /**
   *  Weights for the W to leptons decays.
   */
  vector<double> _wleptonwgt;
  //@}
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of SMWZDecayer.
 */
template <>
 struct BaseClassTrait<Herwig::SMWZDecayer,1> {
    /** Typedef of the base class of SMWZDecayer. */
   typedef Herwig::DecayIntegrator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
 struct ClassTraits<Herwig::SMWZDecayer>
  : public ClassTraitsBase<Herwig::SMWZDecayer> {
   /** Return the class name.*/
  static string className() { return "Herwig::SMWZDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwPerturbativeDecay.so"; }

};

/** @endcond */

}

#include "SMWZDecayer.icc"

#endif /* HERWIG_SMWZDecayer_H */
