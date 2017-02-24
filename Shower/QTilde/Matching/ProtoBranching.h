// -*- C++ -*-
#ifndef HERWIG_ProtoBranching_H
#define HERWIG_ProtoBranching_H
//
// This is the declaration of the ProtoBranching class.
//

namespace Herwig {

using namespace ThePEG;

/**
 *  Declare the pointers
 */
class ProtoBranching;
ThePEG_DECLARE_POINTERS(Herwig::ProtoBranching,ProtoBranchingPtr);

/**
 *  Class to store a prototype branching 
 */
class ProtoBranching : public Base {

public:

  /**
   *  Default constructor
   */
  ProtoBranching() {}

  /**
   *  Constructor
   */
  ProtoBranching(tcPDPtr part, HardBranching::Status status,
		 const Lorentz5Momentum & momentum,
		 tSudakovPtr sudakov)
    : part_(part), status_(status), momentum_(momentum),
      sudakov_(sudakov), type_(ShowerPartnerType::Undefined)
  {}

  /**
   *  Id of the brnaching particle
   */
  long id() { return part_->id();}

  /**
   *  The ParticleData
   */
  tcPDPtr particle() {return part_;}

  /**
   *  Status of the branching
   */
  HardBranching::Status status() {return status_;}

  /**
   *  Set the parent
   */
  tProtoBranchingPtr parent() {return parent_;}

  /**
   *  Get the parent
   */
  void parent(tProtoBranchingPtr in) {parent_=in;}

  /**
   *  Children
   */
  vector<tProtoBranchingPtr> children() {return children_;}

  /**
   *  Add a child
   */
  void addChild(tProtoBranchingPtr in ) {children_.push_back(in);}

  /**
   *  Back children
   */
  vector<tProtoBranchingPtr> backChildren() {return backChildren_;}

  /**
   *  Add a child
   */
  void addBackChild(tProtoBranchingPtr in ) {backChildren_.push_back(in);}

  /**
   *  momentum
   */
  const Lorentz5Momentum & momentum() {return momentum_;}

  /**
   *  Get the Sudakov
   */
  tSudakovPtr sudakov() {return sudakov_;}

  /**
   *  Set the Sudakov
   */
  void sudakov(tSudakovPtr in) { sudakov_=in; }

  /**
   *  Type of branching
   */
  ShowerPartnerType type() const {
    return type_;
  }

  /**
   *  Type of branching
   */
  void type(ShowerPartnerType in) {
    type_ = in;
    assert(type_!=ShowerPartnerType::Undefined);
  }

  /**
   *  Colour line
   */
  tColinePtr colourLine() const {
    return colourLine_;
  }
  /**
   *  Anticolour line
   */
  tColinePtr antiColourLine() const {
    return antiColourLine_;
  }

  /**
   *  Colour line
   */
  void colourLine(tColinePtr in) {
    colourLine_ = in;
  }
  /**
   *  Anticolour line
   */
  void antiColourLine(tColinePtr in) {
    antiColourLine_ = in;
  }

private:

  /**
   *  PDG code
   */
  tcPDPtr part_;

  /**
   *  status
   */
  HardBranching::Status status_;

  /**
   *  Momentum
   */
  Lorentz5Momentum momentum_;

  /**
   *  Sudakov
   */
  tSudakovPtr sudakov_;

  /**
   *  children
   */
  vector<tProtoBranchingPtr> children_;

  /**
   *  back children
   */
  vector<tProtoBranchingPtr> backChildren_;

  /**
   *  parent
   */
  tProtoBranchingPtr parent_;

  /**
   *  The type of branching
   */
  ShowerPartnerType type_;

  /**
   *   Colour lines
   */
  //@{
  /**
   *  Colour line
   */
  tColinePtr colourLine_;

  /**
   *  Anticolour line
   */
  tColinePtr antiColourLine_;
  //@}

};

}

#endif /* HERWIG_ProtoBranching_H */
