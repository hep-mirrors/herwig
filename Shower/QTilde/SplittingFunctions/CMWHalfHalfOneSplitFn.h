  // -*- C++ -*-
  //
  // CMWHalfHalfOneSplitFn.h is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2019 The Herwig Collaboration
  //
  // Herwig is licenced under version 3 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
#ifndef HERWIG_CMWHalfHalfOneSplitFn_H
#define HERWIG_CMWHalfHalfOneSplitFn_H
  //
  // This is the declaration of the CMWHalfHalfOneSplitFn class.
  //

#include "HalfHalfOneSplitFn.h"
#include "Herwig/Shower/ShowerAlpha.h"


namespace Herwig {
  
  using namespace ThePEG;
  
  /** \ingroup Shower
   *
   * This class provides the concrete implementation
   * of the CMW enhanced expressions for the
   * splitting function for \f$\frac12\to q\frac12 1\f$.
   *
   * The kernel uses the same overestimate as the
   * corresponding HalfHalfOneSplitFn and thus only needs to
   * implement the spitting function and ratio to the overestimate.
   *
   * TODO: For a more efficient sampling one needs can rewrite the
   *       overestimation to contain the alpha_max*Kgmax factors.
   *
   * @see \ref CMWHalfHalfOneSplitFnInterfaces "The interfaces"
   * defined for CMWHalfHalfOneSplitFn.
   */
  class CMWHalfHalfOneSplitFn: public HalfHalfOneSplitFn {
    
  public:
    
    /**
     *   Methods to return the splitting function.
     */
      //@{
    /**
     *  Very similar to HalfHalfOneSplitFn.
     *  Here the kernel only contains the soft part multiplied by the 
     *  alphas/2pi * Kg from 
     *  Nucl.Phys. B349 (1991) 635-654
     *
     */
    virtual double P(const double z, const Energy2 t, const IdList & ids,
                     const bool mass, const RhoDMatrix & rho) const;
    
    /**
     *  Very similar to HalfHalfOneSplitFn.
     *  Since we use only the 1/1-z part for overestimating the kernel 
     *  in the first place we can keep the same overestimation related functions
     *  for the CMW kernels.
     */
    virtual double ratioP(const double z, const Energy2 t, const IdList & ids,
                          const bool mass, const RhoDMatrix & rho) const;
    
    /**
     * Return the correction term from:
     * Nucl.Phys. B349 (1991) 635-654
     */
    double Kg(Energy2 )const{
        //TODO: Should be t-dependent
      int Nf=5;//alpha_->Nf(t)
      return (3.*(67./18.-1./6.*sqr(Constants::pi))-5./9.*Nf);
    }
    
    //@}
    
    
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
    virtual void doinit(){
      HalfHalfOneSplitFn::doinit();
    };
      //@}
    
    
  private:
    
      // Pointer to the alpha_s object in use.
    ShowerAlphaPtr alpha_;
      // Provide information if the kernel is used for initial state.
      // as the pt definition contains an additional factor of z.
    bool isIS_=false;
    
  protected:
    
    /** @name Clone Methods. */
      //@{
    /**
     * Make a simple clone of this object.
     * @return a pointer to the new object.
     */
    virtual IBPtr clone() const {return new_ptr(*this);}
    
    /** Make a clone of this object, possibly modifying the cloned object
     * to make it sane.
     * @return a pointer to the new object.
     */
    virtual IBPtr fullclone() const {return new_ptr(*this);}
      //@}
    
  private:
    
    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    CMWHalfHalfOneSplitFn & operator=(const CMWHalfHalfOneSplitFn &) = delete;
    
    
  };
  
}

#endif /* HERWIG_CMWHalfHalfOneSplitFn_H */
