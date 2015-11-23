// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HiggsVBFProcessConstructor class.
//

#include "HiggsVBFProcessConstructor.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/MatrixElement/General/GeneralfftoffH.h"

using namespace Herwig;

HiggsVBFProcessConstructor::HiggsVBFProcessConstructor()
  : _type(true), _shapeOpt(1), _intermediates(0) {
}

IBPtr HiggsVBFProcessConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr HiggsVBFProcessConstructor::fullclone() const {
  return new_ptr(*this);
}

void HiggsVBFProcessConstructor::persistentOutput(PersistentOStream & os) const {
  os << _higgs << _type << _shapeOpt;
}

void HiggsVBFProcessConstructor::persistentInput(PersistentIStream & is, int) {
  is >> _higgs >> _type >> _shapeOpt;
}

ClassDescription<HiggsVBFProcessConstructor> 
HiggsVBFProcessConstructor::initHiggsVBFProcessConstructor;
// Definition of the static class description member.

void HiggsVBFProcessConstructor::Init() {

  static ClassDocumentation<HiggsVBFProcessConstructor> documentation
    ("The HiggsVBFProcessConstructor class generates hard processes for"
     " Higgs boson production in association with a vector boson in general models.");

  static RefVector<HiggsVBFProcessConstructor,ParticleData> interfaceHiggsBoson
    ("HiggsBoson",
     "The possible Higgs bosons",
     &HiggsVBFProcessConstructor::_higgs, -1, false, false, true, false, false);

  static Switch<HiggsVBFProcessConstructor,unsigned int> interfaceShapeOption
    ("ShapeScheme",
     "Option for the treatment of the Higgs resonance shape",
     &HiggsVBFProcessConstructor::_shapeOpt, 2, false, false);
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

  static Switch<HiggsVBFProcessConstructor,unsigned int> interfaceProcess
    ("Process",
     "Which processes to include",
     &HiggsVBFProcessConstructor::_intermediates, 0, false, false);
  static SwitchOption interfaceProcessBoth
    (interfaceProcess,
     "Both",
     "Include both WW and ZZ processes",
     0);
  static SwitchOption interfaceProcessWW
    (interfaceProcess,
     "WW",
     "Only include WW processes",
     1);
  static SwitchOption interfaceProcessZZ
    (interfaceProcess,
     "ZZ",
     "Only include ZZ processes",
     2);

  static Switch<HiggsVBFProcessConstructor,bool> interfaceCollisionType
    ("CollisionType",
     "Type of collision",
     &HiggsVBFProcessConstructor::_type, true, false, false);
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

}

void HiggsVBFProcessConstructor::constructDiagrams() {
  if(_higgs.empty() || !subProcess() ) return;
  tPDPtr Wplus  = getParticleData(ParticleID::Wplus);
  tPDPtr Wminus = getParticleData(ParticleID::Wminus);
  tPDPtr Z0     = getParticleData(ParticleID::Z0);
  for(unsigned int ix=0;ix<_higgs.size();++ix)
    _higgs[ix]->init();
  for(unsigned int ix=0;ix<2;++ix) {
    if( ( ix == 0 && _intermediates == 2 ) ||
	( ix == 1 && _intermediates == 1 )) continue;
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
	// check incoming W+W- or ZZ allowed
	if(ix==0) {
	  if(!vertex->isIncoming(Wminus)||
	     !vertex->isIncoming(Wplus)) continue;
	}
	else {
	  if(!vertex->isIncoming(Z0)) continue;
	}
	// check outgoing higgs allowed
	if(!vertex->isOutgoing(*ih)) continue;
 	// create the MatrixElement object
 	string objectname ("/Herwig/MatrixElements/");
 	string classname("Herwig::GeneralfftoffH");
 	if(_type) objectname += "MEPP2";
 	else      objectname += "MEee2";
	string bos = ix==0 ? "W+W+" : "ZOZO";
	objectname += bos;
 	objectname += (**ih).PDGName();
	GeneralfftoffHPtr matrixElement = dynamic_ptr_cast<GeneralfftoffHPtr>
	  (generator()->preinitCreate(classname, objectname));
	if( !matrixElement )
	  throw Exception()
	    << "HiggsVBFProcessConstructor::constructDiagrams() "
	    << " Failed to construct matrix element for "
	    << bos  << " + "
	    << (**ih).PDGName() << " production"
	    << Exception::runerror;
	GeneralfftoffH::Process process = _type ? 
	  GeneralfftoffH::Hadron : GeneralfftoffH::Lepton;
	// set the information
	matrixElement->setProcessInfo( process, *ih, svert,_shapeOpt,
				       ix+1 );
	// insert it
	generator()->preinitInterface(subProcess(), "MatrixElements", 
				      subProcess()->MEs().size(),
				      "insert", matrixElement->fullName());
      }
    }
  }
}
