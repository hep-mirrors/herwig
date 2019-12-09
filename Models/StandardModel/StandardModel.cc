// -*- C++ -*-
//
// StandardModel.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StandardModel class.
//

#include "StandardModel.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/General/ModelGenerator.h"
#include "ThePEG/Repository/BaseRepository.h"

using namespace Herwig;

StandardModel::StandardModel() {}

StandardModel::~StandardModel() {}

StandardModel::StandardModel(const StandardModel & x)
  : StandardModelBase(x), 
    FFZVertex_ (x.FFZVertex_),
    FFPVertex_ (x.FFPVertex_) , FFGVertex_ (x.FFGVertex_) ,
    FFWVertex_ (x.FFWVertex_) , FFHVertex_ (x.FFHVertex_) ,
    WWHVertex_ (x.WWHVertex_) ,
    GGGVertex_ (x.GGGVertex_) ,
    WWWVertex_ (x.WWWVertex_) , GGGGVertex_(x.GGGGVertex_),
    WWWWVertex_(x.WWWWVertex_), HGGVertex_ (x.HGGVertex_) ,
    HPPVertex_ (x.HPPVertex_) , HHHVertex_ (x.HHHVertex_) ,
    WWHHVertex_ (x.WWHHVertex_) ,
    vertexList_(x.vertexList_), extraVertices_(x.extraVertices_),
    runningMass_(x.runningMass_),modelGenerator_(x.modelGenerator_),
    couplings_(x.couplings_)
{}

IBPtr StandardModel::clone() const {
  return new_ptr(*this);
}

IBPtr StandardModel::fullclone() const {
  return new_ptr(*this);
}

void StandardModel::doinit() {
  if(runningMass_) {
    runningMass_->init();
  }
  //add Standard Model vertices
  if ( registerDefaultVertices() ) {
    addVertex(FFZVertex_);
    addVertex(FFPVertex_);
    addVertex(FFGVertex_);
    addVertex(FFWVertex_);
    addVertex(vertexFFH());
    addVertex(vertexWWH());
    addVertex(GGGVertex_);
    addVertex(WWWVertex_);
    addVertex(GGGGVertex_);
    addVertex(WWWWVertex_);
    addVertex(vertexHGG());
    addVertex(HPPVertex_);
    if(HHHVertex_ ) addVertex(HHHVertex_);
    if(WWHHVertex_) addVertex(WWHHVertex_);
  }
  if(couplings_.find("QED")==couplings_.end()) {
    couplings_["QED"] = make_pair(1,99);
  }
  if(couplings_.find("QCD")==couplings_.end()) {
    couplings_["QCD"] = make_pair(2,99);
  }
  StandardModelBase::doinit();
}

void StandardModel::persistentOutput(PersistentOStream & os) const {
  os << FFZVertex_ <<FFPVertex_ << FFGVertex_ << FFWVertex_ 
     << FFHVertex_ << WWHVertex_ << GGGGVertex_ << WWWWVertex_
     << GGGVertex_ << WWWVertex_  << HGGVertex_  << HPPVertex_ 
     << HHHVertex_ << WWHHVertex_ 
     << runningMass_ << vertexList_ << extraVertices_ << modelGenerator_;
}

void StandardModel::persistentInput(PersistentIStream & is, int) {
  is >> FFZVertex_ >> FFPVertex_ >> FFGVertex_ >> FFWVertex_
     >> FFHVertex_ >> WWHVertex_ >> GGGGVertex_ >> WWWWVertex_
     >> GGGVertex_ >> WWWVertex_ >> HGGVertex_  >> HPPVertex_ 
     >> HHHVertex_ >> WWHHVertex_ 
     >> runningMass_ >> vertexList_ >> extraVertices_ >> modelGenerator_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<StandardModel,StandardModelBase>
describeHerwigStandardModel("Herwig::StandardModel", "Herwig.so");

void StandardModel::Init() {

  static Reference<StandardModel,AbstractFFVVertex> interfaceVertexFFZ
    ("Vertex/FFZ",
     "Reference to the Standard Model FFZ Vertex",
     &StandardModel::FFZVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractFFVVertex> interfaceVertexFFP
    ("Vertex/FFP",
     "Reference to the Standard Model FFP Vertex",
     &StandardModel::FFPVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractFFVVertex> interfaceVertexFFG
    ("Vertex/FFG",
     "Reference to the Standard Model FFG Vertex",
     &StandardModel::FFGVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractFFVVertex> interfaceVertexFFW
    ("Vertex/FFW",
     "Reference to the Standard Model FFW Vertex",
     &StandardModel::FFWVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractFFSVertex> interfaceVertexFFH
    ("Vertex/FFH",
     "Reference to the Standard Model FFH Vertex.",
     &StandardModel::FFHVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractVVVVertex> interfaceVertexGGG
    ("Vertex/GGG",
     "Reference to the Standard Model GGG Vertex",
     &StandardModel::GGGVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractVVVVertex> interfaceVertexWWW
    ("Vertex/WWW",
     "Reference to the Standard Model WWW Vertex",
     &StandardModel::WWWVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractVVSVertex> interfaceVertexWWH
    ("Vertex/WWH",
     "Reference to the Standard Model WWH Vertex",
     &StandardModel::WWHVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractVVVVVertex> interfaceVertexWWWW
    ("Vertex/WWWW",
     "Reference to the Standard Model WWWW Vertex",
     &StandardModel::WWWWVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractVVVVVertex> interfaceVertexGGGG
    ("Vertex/GGGG",
     "Reference to the Standard Model GGGG Vertex",
     &StandardModel::GGGGVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractVVSVertex> interfaceVertexHGG
    ("Vertex/HGG",
     "Reference to the StandardModel HGG Vertex",
     &StandardModel::HGGVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractVVSVertex> interfaceVertexHPP
    ("Vertex/HPP",
     "Reference to StandardModel HPPVertex",
     &StandardModel::HPPVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractSSSVertex> interfaceVertexHHH
    ("Vertex/HHH",
     "Reference to the Standard Model HHHVertex",
     &StandardModel::HHHVertex_, false, false, true, true);

  static Reference<StandardModel,AbstractVVSSVertex> interfaceVertexWWHH
    ("Vertex/WWHH",
     "Reference to the Standard Model WWHHVertex",
     &StandardModel::WWHHVertex_, false, false, true, true);

  static RefVector<StandardModel,VertexBase> interfaceExtraVertices
    ("ExtraVertices",
     "Additional vertices to be considered in automatic ME construction.",
     &StandardModel::extraVertices_, -1, true, false, true, false, false);

  static Reference<StandardModel,RunningMassBase> interfaceRunningMass
    ("RunningMass",
     "Reference to the running mass object",
     &StandardModel::runningMass_, false, false, true, false);
  
  static Reference<StandardModel,Herwig::ModelGenerator> interfaceModelGenerator
    ("ModelGenerator",
     "Pointer to ModelGenerator class",
     &StandardModel::modelGenerator_, false, false, true, true);

  static ClassDocumentation<StandardModel> documentation
    ("The StandardModel class inherits from StandardModelBase"
     "and supplies additional couplings and access to the StandardModel"
     "vertices for helicity amplitude calculations" );

}

void StandardModel::resetMass(long id, Energy mass,tPDPtr part) {
  if(!part) part = getParticleData(id);
  if(!part) return;
  const InterfaceBase * ifb = BaseRepository::FindInterface(part, "NominalMass");
  ostringstream os;
  os << setprecision(12) << abs(mass/GeV);
  ifb->exec(*part, "set", os.str());
}
