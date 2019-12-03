// -*- C++ -*-
//
// ADDModel.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ADDModel class.
//

#include "ADDModel.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

void ADDModel::doinit() {
  addVertex(FFGRVertex_);
  addVertex(VVGRVertex_);
  addVertex(SSGRVertex_);
  addVertex(FFGGRVertex_);
  addVertex(FFWGRVertex_);
  addVertex(GGGGRVertex_);
  addVertex(WWWGRVertex_);
  BSMModel::doinit();
}

void ADDModel::persistentOutput(PersistentOStream & os) const {
  os << ounit(mPlanckBar_,GeV) << ounit(md_,GeV) << delta_
     << ounit(lambdaT_,GeV)
     << FFGRVertex_ << VVGRVertex_ << SSGRVertex_ 
     << FFGGRVertex_ << FFWGRVertex_ 
     << GGGGRVertex_ << WWWGRVertex_;
}

void ADDModel::persistentInput(PersistentIStream & is, int) {
  is >> iunit(mPlanckBar_,GeV) >> iunit(md_,GeV) >> delta_
     >> iunit(lambdaT_,GeV)
     >> FFGRVertex_ >> VVGRVertex_ >> SSGRVertex_
     >> FFGGRVertex_ >> FFWGRVertex_ 
     >> GGGGRVertex_ >> WWWGRVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ADDModel,BSMModel>
describeHerwigADDModel("Herwig::ADDModel", "HwADDModel.so");

void ADDModel::Init() {
  
  static Reference<ADDModel,ThePEG::Helicity::AbstractFFTVertex> interfaceVertexFFGR
    ("Vertex/FFGR",
     "Reference to the fermion-fermion-graviton vertex",
     &ADDModel::FFGRVertex_, false, false, true, false, false);

  static Reference<ADDModel,ThePEG::Helicity::AbstractVVTVertex> interfaceVertexVVGR
    ("Vertex/VVGR",
     "Reference to the vector-vector-graviton vertex",
     &ADDModel::VVGRVertex_, false, false, true, false, false);

  static Reference<ADDModel,ThePEG::Helicity::AbstractSSTVertex> interfaceVertexSSGR
    ("Vertex/SSGR",
     "Reference to the scalar-scalar-graviton vertex",
     &ADDModel::SSGRVertex_, false, false, true, false, false);
  
  static Reference<ADDModel,ThePEG::Helicity::AbstractFFVTVertex> interfaceVertexFFGGR
    ("Vertex/FFGGR",
     "Reference to the fermion-antifermion-gluon graviton vertex",
     &ADDModel::FFGGRVertex_, false, false, true, false, false);
  
  static Reference<ADDModel,ThePEG::Helicity::AbstractFFVTVertex> interfaceVertexFFWGR
    ("Vertex/FFWGR",
     "Reference to the fermion-antifermion-weak vector boson graviton vertex",
     &ADDModel::FFWGRVertex_, false, false, true, false, false);
  
  static Reference<ADDModel,ThePEG::Helicity::AbstractVVVTVertex> interfaceVertexGGGGR
    ("Vertex/GGGGR",
     "Reference to the three gluon graviton vertex",
     &ADDModel::GGGGRVertex_, false, false, true, false, false);
  
  static Reference<ADDModel,ThePEG::Helicity::AbstractVVVTVertex> interfaceVertexWWWGR
    ("Vertex/WWWGR",
     "Reference to the three weak vector boson graviton vertex",
     &ADDModel::WWWGRVertex_, false, false, true, false, false);
  
  static ClassDocumentation<ADDModel> documentation
    ("The ADDModel class replaces the Standard Model class for the"
     " ADD model");
  
  static Parameter<ADDModel,unsigned int> interfaceDelta
    ("Delta",
     "Number of extra dimensions",
     &ADDModel::delta_, 2, 2, 1000,
     false, false, Interface::limited);

  static Parameter<ADDModel,Energy> interfaceReducedPlanckMass
    ("Reduced4dPlanckMass",
     "The reduced planck mass in 4 dimensions",
     &ADDModel::mPlanckBar_, GeV, 2.4e18*GeV, 1e17*GeV, 1e20*GeV,
     false, false, Interface::limited);

  static Parameter<ADDModel,Energy> interfaceDdPlanckMass
    ("DdPlanckMass",
     "The d dimension planck mass",
     &ADDModel::md_, GeV, 1000.*GeV, 100.0*GeV, 1e6*GeV,
     false, false, Interface::limited);

  static Parameter<ADDModel,Energy> interfaceLambdaT
    ("LambdaT",
     "The cut-off for virtual graviton processes",
     &ADDModel::lambdaT_, GeV, 1000.*GeV, 100.*GeV, 100000.0*GeV,
     false, false, Interface::limited);

}
