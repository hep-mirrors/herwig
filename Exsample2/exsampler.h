// -*- C++ -*-

// exsample2/exsampler.h is part of matchbox
// (C) 2008 Simon Platzer -- sp@particle.uni-karlsruhe.de

#ifndef matchbox_exsample2_exsampler_h
#define matchbox_exsample2_exsampler_h

#include <map>

#include "ThePEG/Handlers/SamplerBase.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Utilities/Exception.h"
#include "ThePEG/Repository/UseRandom.h"

#include "exsample/generator.h"

namespace matchbox {

  /**
   * exsample2 function interface conforming
   * wrapper around standard event handler
   */
  struct eg_exsample2_wrapper {

    inline eg_exsample2_wrapper ()
      : bin_(-1), event_handler_() {}

    inline eg_exsample2_wrapper (int b,
				 ThePEG::tStdEHPtr eh)
      : bin_(b), event_handler_(eh) {}
				   

    /**
     * Evaluate with given random numbers.
     */
    inline double evaluate (const std::vector<double>& p) const {
      double ret;
      try {
	ret = event_handler_->dSigDR(p) / ThePEG::nanobarn;
      } catch (ThePEG::Veto&) {
	ret = 0.0;
      } catch (...) {
	throw;
      }
      return ret;
    }

    /**
     * Return the dimensionality in random
     * numbers needed to generate events
     */
    inline std::size_t dimension () const {
      return event_handler_->nDim(bin_);
    }

    /**
     * Return the lower left and upper right
     * corners of the support of this function
     */
    inline std::pair<std::vector<double>,std::vector<double> > support () const {
      std::vector<double> lower(dimension(),0.);
      std::vector<double> upper(dimension(),1.);
      return std::make_pair(lower,upper);
    }

    /**
     * Indicate start of presampling
     */
    inline void start_presampling () { }

    /**
     * Indicate end of presampling
     */
    inline void stop_presampling () { }

  private:

    int bin_;
    ThePEG::tStdEHPtr event_handler_;

  };

  /**
   * Here is the documentation of the exsampler class.
   *
   * @see \ref exsamplerInterfaces "The interfaces"
   * defined for exsampler.
   */
  class exsampler
    : public ThePEG::SamplerBase {

  public:

    /** @name Standard constructors and destructors. */
    //@{
    /**
     * The default constructor.
     */
    exsampler();

    /**
     * The destructor.
     */
    virtual ~exsampler();
    //@}

  public:

    typedef exsample::generator<eg_exsample2_wrapper,ThePEG::UseRandom> sampler_type;

    /** @name Virtual functions from SamplerBase. */
    //@{
    /**
     * Initialize the the sampler, possibly doing presampling of the
     * phase space.
     */
    virtual void initialize();

    /**
     * Generarate a new phase space point and return a weight associated
     * with it. This weight should preferably be 1.
     */
    virtual double generate();

    /**
     * Reject the last chosen phase space point.
     */
    virtual inline void rejectLast() {
      samplers_[last_bin_].reject();
      sum_weights_ -= last_weight_;
    }

    /**
     * If the sampler is able to sample several different functions
     * separately, this function should return the last chosen
     * function. This default version always returns 0.
     */
    virtual inline int lastBin() const { return last_bin_; }

    /**
     * Return the total integrated cross section determined from the
     * Monte Carlo sampling so far.
     */
    virtual inline ThePEG::CrossSection integratedXSec() const {
      return integral_ * ThePEG::nanobarn;
    }

    /**
     * Return the error on the total integrated cross section determined
     * from the Monte Carlo sampling so far.
     */
    virtual ThePEG::CrossSection integratedXSecErr() const {
      return std::sqrt(variance_) * ThePEG::nanobarn;
    }

    /**
     * Return the overestimated integrated cross section.
     */
    virtual ThePEG::CrossSection maxXSec() const {
      return max_integral_ * ThePEG::nanobarn;
    }

    /**
     * Return the sum of the weights returned by generate() so far (of
     * the events that were not rejeted).
     */
    virtual inline double sumWeights() const { return sum_weights_; }
    //@}

  public:

    /** @name Functions used by the persistent I/O system. */
    //@{
    /**
     * Function used to write out object persistently.
     * @name os the persistent output stream written to.
     */
    void persistentOutput(ThePEG::PersistentOStream & os) const;

    /**
     * Function used to read in object persistently.
     * @name is the persistent input stream read from.
     * @name version the version number of the object when written.
     */
    void persistentInput(ThePEG::PersistentIStream & is, int version);
    //@}

    /**
     * The standard Init function used to initialize the interfaces.
     * Called exactly once for each class by the class description system
     * before the main function starts or
     * when this class is dynamically loaded.
     */
    static void Init();

  public:

    struct maxtry_exception : public ThePEG::Exception {};

  protected:

    /** @name Standard Interfaced functions. */
    //@{

    /**
     * Initialize this object. Called in the run phase just before
     * a run begins.
     */
    virtual void doinitrun();

    /**
     * Finalize this object. Called in the run phase just after a
     * run has ended. Used eg. to write out statistics.
     */
    virtual void dofinish();

    //@}

  protected:

    /** @name Clone Methods. */
    //@{
    /**
     * Make a simple clone of this object.
     * @return a pointer to the new object.
     */
    virtual ThePEG::IBPtr clone() const;

    /** Make a clone of this object, possibly modifying the cloned object
     * to make it sane.
     * @return a pointer to the new object.
     */
    virtual ThePEG::IBPtr fullclone() const;
    //@}

  private:

    /// update integrals and stuff
    void update ();

    /// get process by bin id
    std::string process(int bin) const;

    ///@name interface helpers to adaption_info objects
    //@{

    unsigned long presampling_points_;

    unsigned long freeze_grid_;

    double efficiency_threshold_;

    double gain_threshold_;

    //@}

    /// map bin ids to generators
    std::vector<sampler_type> samplers_;

    /// map bin ids to wrappers used
    std::vector<eg_exsample2_wrapper> eg_wrappers_;

    /// map used to select bins
    std::map<double,int> bin_selector_;

    /// current integrals of the bins
    std::vector<double> integrals_;

    /// map bins to number of missing events
    std::vector<long> missing_events_;

    /// bins which do need oversampling
    std::set<int> oversampling_bins_;

    /// the sum of accepted weights
    double sum_weights_;

    /// the total integrated corss section
    double integral_;

    /// sum of absoulte values of integrals
    double integral_abs_;

    /// the ucertainty on the integrated cross section
    double variance_;

    /// the maximum integrated cross section
    double max_integral_;

    /// the last selected bin
    int last_bin_;

    /// the last generated weight
    double last_weight_;

    /// dummy statistics to fix functionality required by rewrite
    exsample::statistics dummy_;

  private:

    /**
     * The static object used to initialize the description of this class.
     * Indicates that this is an abstract class with persistent data.
     */
    static ThePEG::ClassDescription<exsampler> initexsampler;

    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    exsampler & operator=(const exsampler &);

  };

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

  /** @cond TRAITSPECIALIZATIONS */

  /** This template specialization informs ThePEG about the
   *  base classes of exsampler. */
  template <>
  struct BaseClassTrait<matchbox::exsampler,1> {
    /** Typedef of the first base class of exsampler. */
    typedef ThePEG::SamplerBase NthBase;
  };

  /** This template specialization informs ThePEG about the name of
   *  the exsampler class and the shared object where it is defined. */
  template <>
  struct ClassTraits<matchbox::exsampler>
    : public ClassTraitsBase<matchbox::exsampler> {
    /** Return a platform-independent class name */
    static string className() { return "matchbox::exsampler"; }
    /**
     * The name of a file containing the dynamic library where the class
     * exsampler is implemented. It may also include several, space-separated,
     * libraries if the class exsampler depends on other classes (base classes
     * excepted). In this case the listed libraries will be dynamically
     * linked in the order they are specified.
     */
    static string library() { return "HwExsample2.so"; }
  };

  /** @endcond */

}

#endif /* matchbox_exsample2_exsampler_h */
