// -*- C++ -*-

// couplings/alpha_s.h is part of matchbox
// (C) 2008 Simon Platzer -- sp@particle.uni-karlsruhe.de

#ifndef matchbox_couplings_alpha_s_h
#define matchbox_couplings_alpha_s_h

#include <string>

#include <array>

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/StandardModel/AlphaSBase.h"

#include "gsl.h"

namespace matchbox {

  using namespace ThePEG;

  template<class AlphaS>
  struct solve_lambda_below {
      
    typedef AlphaS alpha_s;
      
    inline solve_lambda_below (alpha_s* a,
			       unsigned int n,
			       Energy2 lambda2n,
			       Energy2 mass2)
      : alpha(a), nf_in(n), lambda2_nf_in(lambda2n), threshold(mass2) {}
      
    alpha_s * alpha;
    unsigned int nf_in;
    Energy2 lambda2_nf_in;
    Energy2 threshold;
      
    inline double operator () (double lambda2) {
      return ((*alpha)(threshold,lambda2_nf_in,nf_in) -
	      (*alpha)(threshold,lambda2*MeV2,nf_in-1));
    }
      
  };
    
  template<class AlphaS>
  struct solve_lambda_above {
      
    typedef AlphaS alpha_s;
      
    inline solve_lambda_above (alpha_s * a,
			       unsigned int n,
			       Energy2 lambda2n,
			       Energy2 mass2)
      : alpha(a), nf_in(n), lambda2_nf_in(lambda2n), threshold(mass2) {}
      
    alpha_s * alpha;
    unsigned int nf_in;
    Energy2 lambda2_nf_in;
    Energy2 threshold;
      
    inline double operator () (double lambda2) {
      return ((*alpha)(threshold,lambda2_nf_in,nf_in) -
	      (*alpha)(threshold,lambda2*MeV2,nf_in+1));
    }
      
  };

  template<class AlphaS>
  struct solve_input_lambda {
      
    typedef AlphaS alpha_s;
      
    inline solve_input_lambda (alpha_s * a,
			       unsigned int n,
			       double inalpha,
			       Energy2 inscale)
      : alpha(a), nf_in(n), alpha_in(inalpha), scale_in(inscale) {}
      
    alpha_s * alpha;
    unsigned int nf_in;
    double alpha_in;
    Energy2 scale_in;
      
    inline double operator () (double lambda2) {
      return ((*alpha)(scale_in,lambda2*MeV2,nf_in) - alpha_in);
    }
      
  };

  /**
   * Base class for the strong coupling.
   *
   * @see \ref alpha_sInterfaces "The interfaces"
   * defined for alpha_s.
   */
  class alpha_s
    : public AlphaSBase {

  public:

    /** @name Standard constructors and destructors. */
    //@{
    /**
     * The default constructor.
     */
    alpha_s();
      
    /**
     * The destructor.
     */
    virtual ~alpha_s();
    //@}

  public:

    /** @name Virtual functions as required by AlphaSBase. */
    //@{
    /**
     * The \f$\alpha_S\f$. Return the QCD coupling for a given \a scale
     * using the given standard model object \a sm.
     */
    virtual inline double value(Energy2 scale, const StandardModelBase &) const {
      return operator() (scale);
    }

    /**
     * Return the flavour thresholds used. The returned vector contains
     * (in position <code>i</code>) the scales when the active number of
     * flavours changes from <code>i</code> to <code>i+1</code>.
     */
    virtual inline vector<Energy2> flavourThresholds() const {
      assert(!nfvector.empty());
      return nfvector;
    }

    /**
     * Return the \f$\Lambda_{QCD}\f$ used for different numbers of
     * active flavours.
     */
    virtual inline vector<Energy> LambdaQCDs() const {
      vector<Energy> res;
      for (size_t k = 0; k < 7; ++k)
	res.push_back(sqrt(lambda_squared_[k]));
      return res;
    }
    //@}      

  public:

    /// return alpha_s as function of scale
    inline double operator () (Energy2 scale) const {

      if ( fixed_ )
	return alpha_s_in_;

      assert(matched());
      unsigned int active = active_flavours(scale_factor_*scale);
      return operator () (scale_factor_*scale,lambda_squared_[active],active);

    }

    /// return alpha_s as function of scale, QCD scale
    /// and number of active flavours
    virtual double operator () (Energy2 scale,
				Energy2 lambda2,
				unsigned int nf) const = 0;

    /// match thresholds and write alpha_s
    /// to specified file; arguments are
    /// Q_low/GeV Q_high/GeV n_steps filename
    string check (string args);

  public:

    /// return minimum number of active flavours
    inline unsigned int min_active_flavours () const { return min_active_flavours_; }

    /// set minimum number of active flavours
    inline void min_active_flavours (unsigned int nf) { min_active_flavours_ = nf; }

    /// return maximum number of active flavours
    inline unsigned int max_active_flavours () const { return max_active_flavours_; }

    /// set maximum number of active flavours
    inline void max_active_flavours (unsigned int nf) { max_active_flavours_ = nf; }

    /// return the number of active flavours at the given scale
    inline unsigned int active_flavours (Energy2 scale) const {
      unsigned int active = 0;
      if (scale > 0.*GeV2) {
	while(quark_mass_squared(active) < scale) {
	  if (++active == max_active_flavours_+1)
	    break;
	}
	active -= 1;
      } else {
	active = 0;
      }
      return active;
    }

    /// return the lambda squared for the given number of flavours
    inline Energy2 lambda_squared (unsigned int f) const {
      assert(f < 7);
      return lambda_squared_[f];
    }

    /// return the mass squared for given flavour
    inline Energy2 quark_mass_squared (unsigned int f) const {
      assert(f < 7);
      return quark_masses_squared_[f];
    }

    /// set the mass squared for given flavour
    inline void quark_mass_squared (unsigned int f, Energy2 m2) {
      assert(f < 7);
      quark_masses_squared_[f] = m2;
      matched_ = false;
    }

  public:

    /// perform the threshold matching
    /// given alpha_s value at reference scale
    void match_thresholds ();

    /// return true, if threshold matching has been
    /// performed
    inline bool matched () const { return matched_; }

  protected:

    /** @name Standard Interfaced functions. */
    //@{

    /**
     * Initialize this object after the setup phase before saving an
     * EventGenerator to disk.
     * @throws InitException if object could not be initialized properly.
     */
    virtual inline void doinit() {
      match_thresholds();
      copy(quark_masses_squared_.begin()+1,
           quark_masses_squared_.end(),nfvector.begin());
      AlphaSBase::doinit();
    }

    //@}

    /// return the scale factor
    double scale_factor () const { return scale_factor_; }

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

  private:

    /**
     * The static object used to initialize the description of this class.
     * Indicates that this is an abstract class with persistent data.
     */
    static AbstractClassDescription<alpha_s> initalpha_s;

    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    alpha_s & operator=(const alpha_s &);

  private:

    unsigned int min_active_flavours_;
    unsigned int max_active_flavours_;

    bool matched_;

    double scale_factor_;

    std::array<Energy2,7> quark_masses_squared_;
    std::array<Energy2,7> lambda_squared_;
    vector<Energy2> nfvector=vector<Energy2>(6);
  
    double alpha_s_in_;
    Energy scale_in_;

    pair<Energy2,Energy2> lambda_range_;

    bool fixed_;

  };

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

  /** @cond TRAITSPECIALIZATIONS */

  /** This template specialization informs ThePEG about the
   *  base classes of alpha_s. */
  template <>
  struct BaseClassTrait<matchbox::alpha_s,1> {
    /** Typedef of the first base class of alpha_s. */
    typedef AlphaSBase NthBase;
  };

  /** This template specialization informs ThePEG about the name of
   *  the alpha_s class and the shared object where it is defined. */
  template <>
  struct ClassTraits<matchbox::alpha_s>
    : public ClassTraitsBase<matchbox::alpha_s> {
    /** Return a platform-independent class name */
    static string className() { return "matchbox::alpha_s"; }
    /**
     * The name of a file containing the dynamic library where the class
     * alpha_s is implemented. It may also include several, space-separated,
     * libraries if the class alpha_s depends on other classes (base classes
     * excepted). In this case the listed libraries will be dynamically
     * linked in the order they are specified.
     */
    static string library() { return "HwDipoleShowerAlphaS.so"; }
  };

  /** @endcond */

}

#endif /* matchbox_couplings_alpha_s_h */
