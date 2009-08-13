// -*- C++ -*-
#ifndef HERWIG_HardBranching_H
#define HERWIG_HardBranching_H
//
// This is the declaration of the HardBranching class.
//

#include "ThePEG/Config/ThePEG.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/SudakovFormFactor.h"
#include "HardBranching.fh"
#include "HardTree.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The HardBranching class is designed to contain the information needed for
 * an individual branching in the POWHEG approach
 */
class HardBranching : public Base {

  /**
   *  The HardTree is friend
   */
  friend class HardTree;

public:

  /**
   *  Enum for the status
   */
  enum Status {Outgoing=0,Incoming,Decay};

public:

  /**
   * The default constructor
   * @param particle The particle which is branching
   * @param sudakov  The Sudakov form factor for the branching
   * @param parent   The parent for the branching
   * @param incoming Whether the particle is incoming or outgoing
   */
  HardBranching(ShowerParticlePtr particle, SudakovPtr sudakov,
		tHardBranchingPtr parent,Status status);

  /**
   * Add a child of the branching
   * @param child The child of the branching
   */
  void addChild(HardBranchingPtr child) {_children.push_back(child);} 

  /**
   *  Clear the children
   */
  void clearChildren() { _children.clear(); }

  /**
   *  Add a backward child
   */
  void addBackChild(HardBranchingPtr child) {
    _back_children.push_back(child);
  }

  /**
   *  Clear the backward children
   */
  void clearBackChildren() { _back_children.clear(); }

  /**
   *  Set the momenta of the particles
   */
  void setMomenta(LorentzRotation R, double alpha, Lorentz5Momentum pt,
		  bool setMomentum=true);

  /**
   *  Use the Sudakov to fix the colours
   */
  void fixColours();

  /**
   *  Set and get members for the private member variables
   */
  //@{
  /**
   * Return the branching particle.
   */
  tShowerParticlePtr branchingParticle() const {return _particle;}

  /**
   * Set the branching particle
   */
  void branchingParticle(ShowerParticlePtr in) {_particle=in;}

  /**
   * Get the original momentum
   */
  const Lorentz5Momentum & original() const {return _original;}

  /**
   * Set the original momentum
   */
  void original(const Lorentz5Momentum & in) {_original=in;}

  /**
   *  Get the p reference vector
   */
  const Lorentz5Momentum & pVector() const {return _p;}

  /**
   *  Set the p reference vector
   */
  void pVector(const Lorentz5Momentum & in) {_p=in;}

  /**
   *  Get the n reference vector
   */
  const Lorentz5Momentum & nVector() const {return _n;}

  /**
   *  Set the n reference vector
   */
  void nVector(const Lorentz5Momentum & in) {_n=in;}

  /**
   *  Get the transverse momentum vector
   */
  const Lorentz5Momentum & qPerp() const {return _qt;}

  /**
   *  Set the transverse momentum vector
   */
  void qPerp(const Lorentz5Momentum & in) {_qt=in;}

  /**
   *  Get the momentum the particle should have as the start of a shower
   */
  const Lorentz5Momentum & showerMomentum() const {return _shower;}

  /**
   *  Set the momentum the particle should have as the start of a shower
   */
  void showerMomentum(const Lorentz5Momentum & in ) {_shower=in;}

  /**
   *  Get the transverse momentum
   */
  Energy pT() const {return _pt;}

  /**
   *  Set the transverse momentum
   */
  void pT(Energy in) { _pt=in;}

  /**
   *  Get whether the branching is incoming, outgoing or decay
   */
  Status status() const {return _status;}

  /**
   *  Set whether the branching is incoming, outgoing or decay
   */
  void status(Status in) {_status=in;}

  /**
   *  The parent of the branching
   */
  tHardBranchingPtr parent() const {return _parent;}

  /**
   *  Set the parent of the branching
   */
  void parent(tHardBranchingPtr in) {_parent=in;}

  /**
   *  The Sudakov form-factor
   */
  SudakovPtr sudakov() const {return _sudakov;}

  /**
   *  The Sudakov form-factor
   */
  void sudakov(SudakovPtr in) {_sudakov=in;}

  /**
   *  set the Sudakov for backward branching
   */
  void backSudakov(SudakovPtr in) {_backSudakov=in;}

  /**
   *  The Sudakov for backwards branching
   */
  SudakovPtr backSudakov() const {return _backSudakov;}

  /**
   *  Get the beam particle
   */
  PPtr beam() const {return _beam;}

  /**
   *  Set the beam particle
   */
  void beam(PPtr in) {_beam=in;}

  /**
   * The children
   */
  vector<HardBranchingPtr> & children() {return _children;}

  /**
   *  The children for backward evolution
   */
  vector<HardBranchingPtr> & backChildren() {
    return _back_children;
  }
  //@}

  /**
   *  Information on the Shower variables for the branching
   */
  //@{
  /**
   *  Get the evolution scale
   */
  Energy scale() const {return _scale;}

  /**
   *  The evolution scale
   */
  void scale(Energy in) {_scale=in;}

  /**
   *  The energy fraction
   */
  double z() const {return _z;}

  /**
   *  The energy fraction
   */
  void z(double in) {_z=in;}

  /**
   *  The azimthual angle
   */
  double phi() const {return _phi;}

  /**
   *  The azimthual angle
   */
  void phi(double in) {_phi=in;}
  //@}

  /**
   *  Colour partners
   */
  //@{
  /**
   *  Get the colour partner
   */
  tHardBranchingPtr colourPartner() const {return _partner;}

  /**
   *  The colour partner of the branching
   */
  void colourPartner(tHardBranchingPtr in) {_partner=in;}

  //@}

private:

  /**
   *  The branching particle
   */
  ShowerParticlePtr _particle;

  /**
   *  The rescaled momentum
   */
  Lorentz5Momentum _original;

  /**
   *  The \f$p\f$ reference vector
   */
  Lorentz5Momentum _p;

  /**
   *  The \f$n\f$ reference vector
   */
  Lorentz5Momentum _n;

  /**
   *  The transverse momentum vector
   */
  Lorentz5Momentum _qt;

  /**
   *  The momentum the particle should have as the start of a shower
   */
  Lorentz5Momentum _shower;

  /**
   *  The transverse momentum
   */
  Energy _pt; 

  /**
   *  Whether the branching is incoming, outgoing or a decay
   */
  Status _status;

  /**
   *  Information on the Shower variables for the branching
   */
  //@{
  /**
   *  The evolution scale
   */
  Energy _scale;

  /**
   *  The energy fraction
   */
  double _z;

  /**
   *  The azimthual angle
   */
  double _phi;
  //@}

  /**
   *  The parent of the branching
   */
  tHardBranchingPtr _parent;

  /**
   *  The Sudakov form-factor
   */
  SudakovPtr _sudakov;

  /**
   *  The Sudakov form-factor for backward branching
   */
  SudakovPtr _backSudakov;

  /**
   * The children
   */
  vector<HardBranchingPtr> _children;

  /**
   * The backward shower children
   */
  vector<HardBranchingPtr> _back_children;

  /**
   *  The beam particle
   */
  PPtr _beam;

  /**
   *  The colour partner
   */
  tHardBranchingPtr _partner;

};

}

#endif /* HERWIG_HardBranching_H */
