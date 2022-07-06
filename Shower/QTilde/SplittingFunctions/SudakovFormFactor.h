// -*- C++ -*-
#ifndef Herwig_SudakovFormFactor_H
#define Herwig_SudakovFormFactor_H
//
// This is the declaration of the SudakovFormFactor class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig/Shower/QTilde/ShowerConfig.h"
#include "Herwig/Shower/QTilde/Kinematics/ShowerKinematics.fh"
#include "ThePEG/EventRecord/RhoDMatrix.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "Herwig/Shower/ShowerAlpha.h"

namespace Herwig {

using namespace ThePEG;

/**  \ingroup Shower
 * Enum to define the possible types of colour structure which can occur in
 * the branching.
 */
enum ColourStructure {Undefined=0,
		      TripletTripletOctet  = 1, OctetOctetOctet       = 2,
		      OctetTripletTriplet  = 3, TripletOctetTriplet   = 4,
		      SextetSextetOctet    = 5, TripletTripletSinglet = 6,
		      OctetOctetSinglet    = 7, Epsilon               = 8,
		      OctetSinglet         = 9,
		      ChargedChargedNeutral=-1,
		      ChargedNeutralCharged=-2,
		      NeutralChargedCharged=-3,
		      EW=-4};

/**
 *  A typedef for the BeamParticleData
 */
typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr;

/**
 * Here is the documentation of the SudakovFormFactor class.
 *
 * @see \ref SudakovFormFactorInterfaces "The interfaces"
 * defined for SudakovFormFactor.
 */
class SudakovFormFactor: public Interfaced {

  /**
   *  The SplittingGenerator is a friend to insert the particles in the 
   *  branchings at initialisation
   */
  friend class SplittingGenerator;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  SudakovFormFactor() : interactionType_(ShowerInteraction::UNDEFINED),
                        angularOrdered_(true),
			colourStructure_(Undefined),
			pdfMax_(35.0), pdfFactor_(0)
  {}
  //@}

public:

  /**
   *  Members to generate the scale of the next branching
   */
  //@{
  /**
   * Return the scale of the next time-like branching. If there is no 
   * branching then it returns ZERO.
   * @param startingScale starting scale for the evolution
   * @param ids The PDG codes of the particles in the splitting
   * @param enhance The radiation enhancement factor
   * defined.
   */
  virtual ShoKinPtr generateNextTimeBranching(const Energy startingScale,
					      const IdList &ids,
					      const RhoDMatrix & rho,
					      double enhance, double detuning) = 0;

  /**
   * Return the scale of the next space-like decay branching. If there is no 
   * branching then it returns ZERO.
   * @param startingScale starting scale for the evolution
   * @param stoppingScale stopping scale for the evolution
   * @param minmass The minimum mass allowed for the spake-like particle.
   * @param ids The PDG codes of the particles in the splitting
   * defined.
   * @param enhance The radiation enhancement factor
   */
  virtual ShoKinPtr generateNextDecayBranching(const Energy startingScale,
						 const Energy stoppingScale,
						 const Energy minmass,
						 const IdList &ids,
						 const RhoDMatrix & rho,
						 double enhance,
						 double detuning) = 0;

  /**
   * Return the scale of the next space-like branching. If there is no 
   * branching then it returns ZERO.
   * @param startingScale starting scale for the evolution
   * @param ids The PDG codes of the particles in the splitting
   * @param x The fraction of the beam momentum
   * defined.
   * @param beam The beam particle
   * @param enhance The radiation enhancement factor
   */
   virtual ShoKinPtr generateNextSpaceBranching(const Energy startingScale,
						 const IdList &ids,double x,
						 const RhoDMatrix & rho,
						 double enhance,
						 tcBeamPtr beam,
						 double detuning) = 0;
  //@}

  /**
   *  Azimuthal angle generation
   */
  //@{
  /**
   * Generate the azimuthal angle of the branching for forward evolution
   * @param particle The branching particle
   * @param ids The PDG codes of the particles in the branchings
   * @param The Shower kinematics
   */
  virtual double generatePhiForward(ShowerParticle & particle,const IdList & ids,
				    ShoKinPtr kinematics, const RhoDMatrix & rho)=0;

  /**
   *  Generate the azimuthal angle of the branching for backward evolution
   * @param particle The branching particle
   * @param ids The PDG codes of the particles in the branchings
   * @param The Shower kinematics
   */
  virtual double generatePhiBackward(ShowerParticle & particle,const IdList & ids,
				     ShoKinPtr kinematics, const RhoDMatrix & rho)=0;

  /**
   *  Generate the azimuthal angle of the branching for ISR in decays
   * @param particle The branching particle
   * @param ids The PDG codes of the particles in the branchings
   * @param The Shower kinematics
   */
  virtual double generatePhiDecay(ShowerParticle & particle,const IdList & ids,
				  ShoKinPtr kinematics, const RhoDMatrix & rho)=0;
  //@}
  
public:

  /**
   *  Return the colour structure
   */
  ColourStructure colourStructure() const {return colourStructure_;}

  /**
   *  Method to check the colours are correct
   */
  bool checkColours(const IdList & ids) const;

  /**
   *  Methods to provide public access to the private member variables
   */
  //@{
  /**
   * Return the pointer to the ShowerAlpha object.
   */
  tShowerAlphaPtr alpha() const { return alpha_; }

public:

  /**
   *  Methods to return the interaction type and order for the Sudakov form factor
   */
  //@{
  /**
   *  Purely virtual method which should determine whether this splitting
   *  function can be used for a given set of particles.
   *  @param ids The PDG codes for the particles in the splitting.
   */
  virtual bool accept(const IdList & ids) const = 0;
  //@}


  /**
   *  Method to return the evolution scale given the
   *  transverse momentum, \f$p_T\f$ and \f$z\f$.
   */
  virtual Energy calculateScale(double z, Energy pt, IdList ids,unsigned int iopt) = 0;
  
  /**
   *  Return the type of the interaction
   */
  ShowerInteraction interactionType() const {return interactionType_;}

  /**
   *  Whether or not the interaction is angular ordered
   */
  bool angularOrdered() const {return angularOrdered_;}

protected:

  /**
   *  Set the particles in the splittings
   */
  void addSplitting(const IdList &);

  /**
   *  Delete the particles in the splittings
   */
  void removeSplitting(const IdList &);

  /**
   *  Access the potential branchings
   */
  const vector<IdList> & particles() const { return particles_;}

  /**
   *  The PDF factor
   */
  unsigned pdfFactor() const {return pdfFactor_;}

  /**
   * Maximum value of the PDF weight
   */
  double pdfMax() const {return pdfMax_;}

  /**
   * Veto on the PDF for the initial-state shower
   * @param t The scale
   * @param x The fraction of the beam momentum
   * @param parton0 Pointer to the particleData for the 
   *                new parent (this is the particle we evolved back to)
   * @param parton1 Pointer to the particleData for the 
   *                original particle
   * @param beam The BeamParticleData object
   */
  bool PDFVeto(const Energy2 t, const double x, const double z,
	       const tcPDPtr parton0, const tcPDPtr parton1,
	       tcBeamPtr beam) const;
  /**
   * The PDF veto ratio
   */
  double PDFVetoRatio(const Energy2 t, const double x, const double z,
               const tcPDPtr parton0, const tcPDPtr parton1,
               tcBeamPtr beam,double factor) const;

  /**
   *   Set the PDF
   */
  void setPDF(tcPDFPtr pdf, Energy scale) {
    pdf_ = pdf;
    freeze_ = scale;
  }

  /**
   *  The veto on the coupling constant
   * @param pt2 The value of ther transverse momentum squared, \f$p_T^2\f$.
   * @return true if vetoed
   */
  bool alphaSVeto(Energy2 pt2) const;

  /**
   * The alpha S veto ratio
   */
  virtual double alphaSVetoRatio(Energy2 pt2,double factor) const;

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SudakovFormFactor & operator=(const SudakovFormFactor &) = delete;

private:

  /**
   *  The interaction type for the splitting function.
   */
  ShowerInteraction interactionType_;

  /**
   *  Whether or not this interaction is angular-ordered
   */
  bool angularOrdered_;

  /**
   *  The colour structure
   */
  ColourStructure colourStructure_;

  /**
   * List of the particles this Sudakov is used for to aid in setting up
   * interpolation tables if needed
   */
  vector<IdList> particles_;

  /**
   *  Stuff for the PDFs
   */
  //@{
  /**
   *  PDf
   */
  tcPDFPtr pdf_;

  /**
   *  Freezing scale
   */
  Energy freeze_;

  /**
   * Maximum value of the PDF weight
   */
  double pdfMax_;

  /**
   *  Option for the inclusion of a factor \f$1/(1-z)\f$ in the PDF estimate
   */
  unsigned pdfFactor_;
  //@}

  /**
   *  Pointer to the coupling for this Sudakov form factor
   */
  ShowerAlphaPtr alpha_;

};

}

#endif /* Herwig_SudakovFormFactor_H */
