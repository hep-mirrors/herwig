// -*- C++ -*-

// couplings/nlo_alpha_s.h is part of matchbox
// (C) 2008 Simon Platzer -- sp@particle.uni-karlsruhe.de

#ifndef matchbox_couplings_nlo_alpha_s_h
#define matchbox_couplings_nlo_alpha_s_h

#include "alpha_s.h"

namespace matchbox {

  using namespace ThePEG;

  /**
   * NLO running alpha_s
   *
   * @see \ref nlo_alpha_sInterfaces "The interfaces"
   * defined for nlo_alpha_s.
   */
  class nlo_alpha_s
    : public alpha_s {

  public:

    /** @name Standard constructors and destructors. */
    //@{
    /**
     * The default constructor.
     */
    nlo_alpha_s();

    /**
     * The destructor.
     */
    virtual ~nlo_alpha_s();
    //@}

  public:

    /// return alpha_s as function of scale, QCD scale
    /// and number of active flavours
    virtual double operator () (Energy2 scale,
				Energy2 lambda2,
				unsigned int nf) const;

    /// return the number of loops which determine this running
    virtual unsigned int nloops () const { return 2; }

  public:

    /** @name Functions used by the persistent I/O system. */
    //@{
    /**
     * Function used to write out object persistently.
     * @name os the persistent output stream written to.
     */
    void persistentOutput(PersistentOStream & os) const;

    /**
     * Function used to read in object persistently.
     * @name is the persistent input stream read from.
     * @name version the version number of the object when written.
     */
    void persistentInput(PersistentIStream & is, int version);
    //@}

    /**
     * The standard Init function used to initialize the interfaces.
     * Called exactly once for each class by the class description system
     * before the main function starts or
     * when this class is dynamically loaded.
     */
    static void Init();

  protected:

    /** @name Standard Interfaced functions. */
    //@{

    /**
     * Initialize this object after the setup phase before saving an
     * EventGenerator to disk.
     * @throws InitException if object could not be initialized properly.
     */
    virtual inline void doinit() {
      freezing_scale_ *= scale_factor();
      alpha_s::doinit();
    }

    //@}

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
     * The static object used to initialize the description of this class.
     * Indicates that this is an abstract class with persistent data.
     */
    static ClassDescription<nlo_alpha_s> initnlo_alpha_s;

    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    nlo_alpha_s & operator=(const nlo_alpha_s &);

  private:

    struct rg_solution {

      inline double operator () (double alpha) {

	double beta0 = (33.-2.*nf)/(12.*Constants::pi);
	double beta1 = (153.-19.*nf)/(24.*sqr(Constants::pi));

	return ((1./alpha)+(beta1/beta0)*log(alpha/(beta0+beta1*alpha))- beta0*slog);

      }

      double slog;
      unsigned int nf;

    };

    Energy freezing_scale_;

    bool exact_evaluation_;

    static rg_solution& rg () {
      static rg_solution rg_;
      return rg_;
    }

    static gsl::bisection_root_solver<rg_solution,100>& rg_solver () {
      static gsl::bisection_root_solver<rg_solution,100> rg_solver_(rg());
      return rg_solver_;
    }

    bool two_largeq_terms_;

  };

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

  /** @cond TRAITSPECIALIZATIONS */

  /** This template specialization informs ThePEG about the
   *  base classes of nlo_alpha_s. */
  template <>
  struct BaseClassTrait<matchbox::nlo_alpha_s,1> {
    /** Typedef of the first base class of nlo_alpha_s. */
    typedef matchbox::alpha_s NthBase;
  };

  /** This template specialization informs ThePEG about the name of
   *  the nlo_alpha_s class and the shared object where it is defined. */
  template <>
  struct ClassTraits<matchbox::nlo_alpha_s>
    : public ClassTraitsBase<matchbox::nlo_alpha_s> {
    /** Return a platform-independent class name */
    static string className() { return "matchbox::nlo_alpha_s"; }
    /**
     * The name of a file containing the dynamic library where the class
     * nlo_alpha_s is implemented. It may also include several, space-separated,
     * libraries if the class nlo_alpha_s depends on other classes (base classes
     * excepted). In this case the listed libraries will be dynamically
     * linked in the order they are specified.
     */
    static string library() { return "HwDipoleShowerAlphaS.so"; }
  };

  /** @endcond */

}

#endif /* matchbox_couplings_nlo_alpha_s_h */
