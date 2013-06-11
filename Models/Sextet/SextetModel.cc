// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SextetModel class.
//

#include "SextetModel.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/Throw.h"

using namespace Herwig;

IBPtr SextetModel::clone() const {
  return new_ptr(*this);
}

IBPtr SextetModel::fullclone() const {
  return new_ptr(*this);
}

void SextetModel::persistentOutput(PersistentOStream & os) const {
  os << VVVVertex_ << VVVVVertex_ << VSSVertex_ << VVSSVertex_
     << FFVVertex_ << FFSVertex_
     << g1L_ << g1R_ << g1pR_ << g1ppR_ << g2_ << g2p_ << g3L_
     << enableScalarSingletY43_ << enableScalarSingletY13_ 
     << enableScalarSingletY23_ << enableScalarTripletY13_ 
     << enableVectorDoubletY16_ << enableVectorDoubletY56_;
}

void SextetModel::persistentInput(PersistentIStream & is, int) {
  is >> VVVVertex_ >> VVVVVertex_ >> VSSVertex_ >> VVSSVertex_
     >> FFVVertex_ >> FFSVertex_
     >> g1L_ >> g1R_ >> g1pR_ >> g1ppR_ >> g2_ >> g2p_ >> g3L_
     >> enableScalarSingletY43_ >> enableScalarSingletY13_ 
     >> enableScalarSingletY23_ >> enableScalarTripletY13_ 
     >> enableVectorDoubletY16_ >> enableVectorDoubletY56_;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<SextetModel,StandardModel>
  describeSextetModel("Herwig::SextetModel", "HwSextetModel.so");

void SextetModel::Init() {

  static ClassDocumentation<SextetModel> documentation
    ("The SextetModel class provides the Model class for models with new scalars"
     " or vectors in the sextet representation of SU(3)");

  static Reference<SextetModel,ThePEG::Helicity::AbstractVVVVertex>
    interfaceVertexVDQVDQG
    ("Vertex/VDQVDQG",
     "The coupling of the gluon to two vector diquarks",
     &SextetModel::VVVVertex_, false, false, true, false, false);

  static Reference<SextetModel,ThePEG::Helicity::AbstractVVVVVertex>
    interfaceVertexVDQVDQGG
    ("Vertex/VDQVDQGG",
     "The coupling of two gluons to two vector diquarks",
     &SextetModel::VVVVVertex_, false, false, true, false, false);

  static Reference<SextetModel,ThePEG::Helicity::AbstractVSSVertex>
    interfaceVertexSDQSDQG
    ("Vertex/SDQSDQG",
     "The coupling of the gluon to two scalar diquarks",
     &SextetModel::VSSVertex_, false, false, true, false, false);

  static Reference<SextetModel,ThePEG::Helicity::AbstractVVSSVertex>
    interfaceVertexSDQSDQGG
    ("Vertex/SDQSDQGG",
     "The coupling of two gluons to two scalar diquarks",
     &SextetModel::VVSSVertex_, false, false, true, false, false);

  static Reference<SextetModel,ThePEG::Helicity::AbstractFFSVertex> 
    interfaceVertexFFSDQ
    ("Vertex/FFSDQ",
     "The coupling of two quarks to the scalar diquark",
     &SextetModel::FFSVertex_, false, false, true, false, false);

  static Reference<SextetModel,ThePEG::Helicity::AbstractFFVVertex> 
    interfaceVertexFFVDQ
    ("Vertex/FFVDQ",
     "The coupling of two quarks to the vector diquark",
     &SextetModel::FFVVertex_, false, false, true, false, false);

  static ParVector<SextetModel,double> interfaceg1L
    ("g1L",
     "The \\f$SU(2)\\f$ quark-doublet coupling to \\f$\\Phi_{6,1,1/3}\\f$.",
     &SextetModel::g1L_, 3, 0.0, 0, 0,
     false, false, Interface::nolimits);

  static ParVector<SextetModel,double> interfaceg1R
    ("g1R",
     "The \\f$SU(2)\\f$ singlet coupling to \\f$\\Phi_{6,1,1/3}\\f$.",
     &SextetModel::g1R_, 3, 0.0, 0, 0,
     false, false, Interface::nolimits);

  static ParVector<SextetModel,double> interfaceg1RPrime
    ("g1RPrime",
     "The \\f$SU(2)\\f$ singlet coupling to \\f$\\Phi_{6,1,1/3}\\f$.",
     &SextetModel::g1pR_, 3, 0.0, 0, 0,
     false, false, Interface::nolimits);

  static ParVector<SextetModel,double> interfaceg1RDoublePrime
    ("g1RDoublePrime",
     "The \\f$SU(2)\\f$ singlet coupling to \\f$\\Phi_{6,1,1/3}\\f$.",
     &SextetModel::g1ppR_, 3, 0.0, 0, 0,
     false, false, Interface::nolimits);

  static ParVector<SextetModel,double> interfaceg2
    ("g2",
     "The coupling to \\f$V^\\mu_{6,2,-1/6}\\f$.",
     &SextetModel::g2_, 3, 0.0, 0, 0,
     false, false, Interface::nolimits);

  static ParVector<SextetModel,double> interfaceg2Prime
    ("g2Prime",
     "The coupling to \\f$V^\\mu_{6,2,5/6}\\f$.",
     &SextetModel::g2p_, 3, 0.0, 0, 0,
     false, false, Interface::nolimits);

  static ParVector<SextetModel,double> interfaceg3L
    ("g3L",
     "Coupling to \\f$\\Phi_{6,3,1/3}\\f$.",
     &SextetModel::g3L_, 3, 0.0, 0, 0,
     false, false, Interface::nolimits);

  static Command<SextetModel> interfaceEnableParticles
    ("EnableParticles",
     "Enable specfic diquarks",
     &SextetModel::doEnable, false);

}

void SextetModel::doinit() {
  StandardModel::doinit();
  if ( !(enableScalarSingletY43_ || enableScalarSingletY13_ 
	 || enableScalarSingletY23_ || enableScalarTripletY13_ 
	 || enableVectorDoubletY16_ || enableVectorDoubletY56_ )) {
    Throw<Exception>() << "You have not enabled any Sextet diquarks. Use e.g.\n"
		       << "  do Model:EnableParticles Scalar Triplet Y=1/3\n"
		       << "to specify the spin, weak isospin and weak hypercharge." 
		       << Exception::runerror;
  }
  addVertex(VVVVertex_);
  addVertex(VVVVVertex_);
  addVertex(VSSVertex_);
  addVertex(VVSSVertex_);
  addVertex(FFVVertex_);
  addVertex(FFSVertex_);
}

string SextetModel::doEnable(string args) {
  int spin=-1;
  int weak=-1;
  int Y[2]={-1000000,-1000000};
  string orig=args;
  while ( !args.empty() ) {
    string arg = StringUtils::car(args);
    args = StringUtils::cdr(args);
    if      ( arg == "Scalar" ) spin=1;
    else if ( arg == "Vector" ) spin=3;
    else if ( arg == "Singlet" ) weak=1;
    else if ( arg == "Doublet" ) weak=2;
    else if ( arg == "Triplet" ) weak=3;
    else {
      if(arg.find("Y=")==string::npos) continue;
      arg = StringUtils::cdr(arg,"=");
      vector<string> split = StringUtils::split(arg,"/");
      if(split.size()!=2) continue;
      istringstream is1(split[0]);
      is1 >> Y[0];
      istringstream is2(split[1]);
      is2 >> Y[1];
    }
  }
  // check we read a value for all three quantum numbers
  if ( spin <0 || weak<0 || 0 || Y[0]== -1000000) {
    return string("SextetModel:EnableParticles couldn't termine spin, weak") + 
      string(" isospin or hypercharge for ") + orig + ".";
  }
  // check the values of Y
  if(!(Y[1]==3||Y[1]==6)) {
    return string("SextetModel:EnableParticles invalid weak") + 
      string(" hypercharge for ") + orig + ".";
  }
  // the various allowed combinations
  bool found = false;
  if(spin == 1 ) {
    found = true;
    if     ( weak == 1 && Y[0] ==  4 && Y[1] == 3) {
      enableScalarSingletY43_ = true;
    }
    else if( weak == 1 && Y[0] ==  1 && Y[1] == 3) {
      enableScalarSingletY13_ = true;
    }
    else if( weak == 1 && Y[0] == -2 && Y[1] == 3) {
      enableScalarSingletY23_ = true;
    }
    else if( weak == 3 && Y[0] ==  1 && Y[1] == 3) {
      enableScalarTripletY13_ = true;
    }
    else
      found = false;
  }
  else if(spin == 3 && weak == 2) {
    found = true;
    if     ( Y[0] == -1 && Y[1] == 6) {
      enableVectorDoubletY16_ = true;
    }
    else if( Y[0] ==  5 && Y[1] == 6) {
      enableVectorDoubletY56_ = true;
    }
    else
      found = false;
  }
  if(!found)
    return string("SextetModel:EnableParticles invalid combination") + 
      string(" of spin, weak isospin or hypercharge for ") + orig + ".";
  else
    return "";
}
