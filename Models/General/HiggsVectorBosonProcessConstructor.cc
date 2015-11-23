// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HiggsVectorBosonProcessConstructor class.
//

#include "HiggsVectorBosonProcessConstructor.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/MatrixElement/General/GeneralfftoVH.h"

using namespace Herwig;

HiggsVectorBosonProcessConstructor::HiggsVectorBosonProcessConstructor()
  : _type(true), _shapeOpt(1) {
}

IBPtr HiggsVectorBosonProcessConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr HiggsVectorBosonProcessConstructor::fullclone() const {
  return new_ptr(*this);
}

void HiggsVectorBosonProcessConstructor::persistentOutput(PersistentOStream & os) const {
  os << _vector << _higgs << _type << _shapeOpt << _alpha;
}

void HiggsVectorBosonProcessConstructor::persistentInput(PersistentIStream & is, int) {
  is >> _vector >> _higgs >> _type >> _shapeOpt >> _alpha;
}

ClassDescription<HiggsVectorBosonProcessConstructor> 
HiggsVectorBosonProcessConstructor::initHiggsVectorBosonProcessConstructor;
// Definition of the static class description member.

void HiggsVectorBosonProcessConstructor::Init() {

  static ClassDocumentation<HiggsVectorBosonProcessConstructor> documentation
    ("The HiggsVectorBosonProcessConstructor class generates hard process for"
     " Higgs boson production in assoication with a vector boson in general models.");

  static RefVector<HiggsVectorBosonProcessConstructor,ParticleData> interfaceVectorBoson
    ("VectorBoson",
     "The possible outgoing vector bosons, must be W/Z",
     &HiggsVectorBosonProcessConstructor::_vector, -1, false, false, true, false, false);

  static RefVector<HiggsVectorBosonProcessConstructor,ParticleData> interfaceHiggsBoson
    ("HiggsBoson",
     "The possible Higgs bosons",
     &HiggsVectorBosonProcessConstructor::_higgs, -1, false, false, true, false, false);

  static Switch<HiggsVectorBosonProcessConstructor,unsigned int> interfaceShapeOption
    ("ShapeScheme",
     "Option for the treatment of the Higgs resonance shape",
     &HiggsVectorBosonProcessConstructor::_shapeOpt, 2, false, false);
  static SwitchOption interfaceStandardShapeFixed
    (interfaceShapeOption,
     "FixedBreitWigner",
     "Breit-Wigner s-channel resonanse",
     1);
  static SwitchOption interfaceStandardShapeRunning
    (interfaceShapeOption,
     "MassGenerator",
     "Use the mass generator to give the shape",
     2);
  static SwitchOption interfaceStandardShapeOn
    (interfaceShapeOption,
     "OnShell",
     "Produce the Higgs on-shell",
     0);

  static Switch<HiggsVectorBosonProcessConstructor,bool> interfaceCollisionType
    ("CollisionType",
     "Type of collision",
     &HiggsVectorBosonProcessConstructor::_type, true, false, false);
  static SwitchOption interfaceCollisionTypeLepton
    (interfaceCollisionType,
     "Lepton",
     "Lepton-Lepton collisions",
     false);
  static SwitchOption interfaceCollisionTypeHadron
    (interfaceCollisionType,
     "Hadron",
     "Hadron-Hadron collisions",
     true);

  static Reference<HiggsVectorBosonProcessConstructor,ShowerAlpha> interfaceAlphaQCD
    ("AlphaQCD",
     "The strong coupling used in the shower for MME or POWHEG corrections.",
     &HiggsVectorBosonProcessConstructor::_alpha, false, false, true, false, false);

}

void HiggsVectorBosonProcessConstructor::constructDiagrams() {
  if(_vector.empty()||_higgs.empty() || !subProcess() ) return;
  // initialise the particles
  for(unsigned int ix=0;ix<_vector.size();++ix)
    _vector[ix]->init();
  for(unsigned int ix=0;ix<_higgs.size();++ix)
    _higgs[ix]->init();
  for(PDVector::const_iterator iv=_vector.begin();
      iv!=_vector.end();++iv) {
    // skip if combination not possible
    if(_type==false && (**iv).id()!=ParticleID::Z0)
      continue;
    else if(_type==true && (abs((**iv).id()) != ParticleID::Wplus &&
			    (**iv).id()      != ParticleID::Z0)) 
      continue;
    // loop over the possible Higgs bosons
    for(PDVector::const_iterator ih=_higgs.begin();
	ih!=_higgs.end();++ih) {
      // check higgs is neutral and scalar
      if((**ih).iCharge()!=0 || (**ih).coloured() ||
	 (**ih).iSpin()!=PDT::Spin0) continue;
      // find a suitable vertex
      for(unsigned int nv = 0; nv < model()->numberOfVertices(); ++nv ) {
	VertexBasePtr vertex = model()->vertex(nv);
	AbstractVVSVertexPtr svert = 
	  dynamic_ptr_cast<AbstractVVSVertexPtr>(vertex);
	if(!svert) continue;
	if(vertex->getNpoint() != 3) continue;
	if(!vertex->isIncoming(*iv)) continue;
	if(!vertex->isOutgoing(*iv)) continue;
	if(!vertex->isOutgoing(*ih)) continue;
	// create the MatrixElement object
	string objectname ("/Herwig/MatrixElements/");
	string classname("Herwig::GeneralfftoVH");
	if(_type) objectname += "MEPP2";
	else      objectname += "MEee2";
	objectname += (**iv).PDGName();
	objectname += (**ih).PDGName();
	GeneralfftoVHPtr matrixElement = dynamic_ptr_cast<GeneralfftoVHPtr>
	  (generator()->preinitCreate(classname, objectname));
	if( !matrixElement )
	  throw Exception()
	    << "HiggsVectorBosonProcessConstructor::constructDiagrams() "
	    << " Failed to construct matrix element for "
	    << (**iv).PDGName() << " + "
	    << (**ih).PDGName() << " production"
	    << Exception::runerror;
	GeneralfftoVH::Process process = GeneralfftoVH::Lepton;
	if(_type) {
	  if((**iv).id()==ParticleID::Z0)
	    process = GeneralfftoVH::HadronZ;
	  else if((**iv).id()==ParticleID::Wplus)
	    process = GeneralfftoVH::HadronWplus;
	  else if((**iv).id()==ParticleID::Wminus)
	    process = GeneralfftoVH::HadronWminus;
	}
	// set the coupling
	generator()->preinitInterface(matrixElement, "Coupling", 
				      "set", _alpha->fullName()); 
	// set the information
	matrixElement->setProcessInfo( process, *ih, svert,_shapeOpt);
	// insert it
	generator()->preinitInterface(subProcess(), "MatrixElements", 
				      subProcess()->MEs().size(),
				      "insert", matrixElement->fullName());
      }
    }
  }
}
