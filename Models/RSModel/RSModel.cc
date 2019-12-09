// -*- C++ -*-
//
// RSModel.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModel class.
//

#include "RSModel.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

void RSModel::doinit() {
  addVertex(FFGRVertex_);
  addVertex(VVGRVertex_);
  addVertex(SSGRVertex_);
  addVertex(FFGGRVertex_);
  addVertex(FFWGRVertex_);
  addVertex(GGGGRVertex_);
  addVertex(WWWGRVertex_);
  BSMModel::doinit();
}

void RSModel::persistentOutput(PersistentOStream & os) const {
  os << ounit(Lambda_pi_,GeV) 
     << FFGRVertex_ << VVGRVertex_ << SSGRVertex_ 
     << FFGGRVertex_ << FFWGRVertex_ 
     << GGGGRVertex_ << WWWGRVertex_;
}

void RSModel::persistentInput(PersistentIStream & is, int) {
  is >> iunit(Lambda_pi_,GeV) 
     >> FFGRVertex_ >> VVGRVertex_ >> SSGRVertex_
     >> FFGGRVertex_ >> FFWGRVertex_ 
     >> GGGGRVertex_ >> WWWGRVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RSModel,BSMModel>
describeHerwigRSModel("Herwig::RSModel", "HwRSModel.so");

void RSModel::Init() {
  

  static Reference<RSModel,ThePEG::Helicity::AbstractFFTVertex> interfaceVertexFFGR
    ("Vertex/FFGR",
     "Reference to the fermion-fermion-graviton vertex",
     &RSModel::FFGRVertex_, false, false, true, false, false);
  
  static Reference<RSModel,ThePEG::Helicity::AbstractVVTVertex> interfaceVertexVVGR
    ("Vertex/VVGR",
     "Reference to the vector-vector-graviton vertex",
     &RSModel::VVGRVertex_, false, false, true, false, false);
  
  static Reference<RSModel,ThePEG::Helicity::AbstractSSTVertex> interfaceVertexSSGR
    ("Vertex/SSGR",
     "Reference to the scalar-scalar-graviton vertex",
     &RSModel::SSGRVertex_, false, false, true, false, false);
  
  static Reference<RSModel,ThePEG::Helicity::AbstractFFVTVertex> interfaceVertexFFGGR
    ("Vertex/FFGGR",
     "Reference to the fermion-antifermion-gluon graviton vertex",
     &RSModel::FFGGRVertex_, false, false, true, false, false);
  
  static Reference<RSModel,ThePEG::Helicity::AbstractFFVTVertex> interfaceVertexFFWGR
    ("Vertex/FFWGR",
     "Reference to the fermion-antifermion-weak vector boson graviton vertex",
     &RSModel::FFWGRVertex_, false, false, true, false, false);
  
  static Reference<RSModel,ThePEG::Helicity::AbstractVVVTVertex> interfaceVertexGGGGR
    ("Vertex/GGGGR",
     "Reference to the three gluon graviton vertex",
     &RSModel::GGGGRVertex_, false, false, true, false, false);
  
  static Reference<RSModel,ThePEG::Helicity::AbstractVVVTVertex> interfaceVertexWWWGR
    ("Vertex/WWWGR",
     "Reference to the three weak vector boson graviton vertex",
     &RSModel::WWWGRVertex_, false, false, true, false, false);
  
  static Parameter<RSModel,Energy> interfaceLambda_pi
    ("Lambda_pi",
     "The coupling of the graviton to matter",
     &RSModel::Lambda_pi_, GeV, 10000*GeV, ZERO, 1.0e12*GeV,
     false, false, false);
  
  static ClassDocumentation<RSModel> documentation
    ("The RSModel class replaces the Standard Model class for the"
     " RS model",
     "The Randall-Sundrum model was constructed from \\cite{Randall:1999ee}.",
     "%\\cite{Randall:1999ee}\n"
     "\\bibitem{Randall:1999ee}\n"
     "  L.~Randall and R.~Sundrum,\n"
     "  ``A large mass hierarchy from a small extra dimension,''\n"
     "  Phys.\\ Rev.\\ Lett.\\  {\\bf 83}, 3370 (1999)\n"
     "  [arXiv:hep-ph/9905221].\n"
     "  %%CITATION = PRLTA,83,3370;%%\n"
     );
  
}
