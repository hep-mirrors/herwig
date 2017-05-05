// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QQHiggsProcessConstructor class.
//

#include "QQHiggsProcessConstructor.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/MatrixElement/General/GeneralQQHiggs.h"

using namespace Herwig;

QQHiggsProcessConstructor::QQHiggsProcessConstructor() 
  : _process(0), _quarkFlavour(0), _shapeOpt(1)
{}

IBPtr QQHiggsProcessConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr QQHiggsProcessConstructor::fullclone() const {
  return new_ptr(*this);
}

void QQHiggsProcessConstructor::persistentOutput(PersistentOStream & os) const {
  os << _process << _quarkFlavour << _higgs << _shapeOpt;
}

void QQHiggsProcessConstructor::persistentInput(PersistentIStream & is, int) {
  is >> _process >> _quarkFlavour >> _higgs >> _shapeOpt;
}

ClassDescription<QQHiggsProcessConstructor> 
QQHiggsProcessConstructor::initQQHiggsProcessConstructor;
// Definition of the static class description member.

void QQHiggsProcessConstructor::Init() {

  static ClassDocumentation<QQHiggsProcessConstructor> documentation
    ("The QQHiggsProcessConstructor class generates hard processes for the"
     " production of the Higgs boson in association with a heavy quark-antiquark"
     " pair in general models.");

  static RefVector<QQHiggsProcessConstructor,ParticleData> interfaceHiggsBoson
    ("HiggsBoson",
     "The possible Higgs bosons",
     &QQHiggsProcessConstructor::_higgs, -1, false, false, true, false, false);

  static Switch<QQHiggsProcessConstructor,unsigned int> interfaceShapeOption
    ("ShapeScheme",
     "Option for the treatment of the Higgs resonance shape",
     &QQHiggsProcessConstructor::_shapeOpt, 2, false, false);
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
  static SwitchOption interfaceStandardShapeYes
    (interfaceShapeOption,
     "OnShell",
     "Produce the Higgs on-shell",
     0);

  static Switch<QQHiggsProcessConstructor,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &QQHiggsProcessConstructor::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all subprocesses",
     0);
  static SwitchOption interfaceProcess1
    (interfaceProcess,
     "gg",
     "Include only gg -> QQbarH processes",
     1);
  static SwitchOption interfaceProcessqbarqbarqbarqbar
    (interfaceProcess,
     "qqbar",
     "Include only qbar qbar -> QQbarH processes",
     2);

  static Switch<QQHiggsProcessConstructor,unsigned int> interfaceQuarkType
    ("QuarkType",
     "The type of quark",
     &QQHiggsProcessConstructor::_quarkFlavour, 6, false, false);
  static SwitchOption interfaceQuarkTypeBottom
    (interfaceQuarkType,
     "Bottom",
     "Produce bottom-antibottom",
     5);
  static SwitchOption interfaceQuarkTypeTop
    (interfaceQuarkType,
     "Top",
     "Produce top-antitop",
     6);
  static SwitchOption interfaceQuarkTypeBottomTop
    (interfaceQuarkType,
     "BottomandTop",
     "Produce bottom-antibottom and top-antitop",
     0);

}


void QQHiggsProcessConstructor::constructDiagrams() {
  if(_higgs.empty() || !subProcess() ) return;
  // initialize the Higgs bosons
  for(unsigned int ix=0;ix<_higgs.size();++ix)
    _higgs[ix]->init();
  long qmin = _quarkFlavour == 0 ? 5 : _quarkFlavour;
  long qmax = _quarkFlavour == 0 ? 6 : _quarkFlavour;
  for(long iq=qmin;iq<=qmax;++iq) {
    tPDPtr qk = getParticleData(iq);
    tPDPtr qb = qk->CC();
    // loop over the possible Higgs bosons
    for(PDVector::const_iterator ih=_higgs.begin();
	ih!=_higgs.end();++ih) {
      // check higgs is neutral and scalar
      if((**ih).iCharge()!=0 || (**ih).coloured() ||
	 (**ih).iSpin()!=PDT::Spin0) continue;
      // find a suitable vertex
      for(unsigned int nv = 0; nv < model()->numberOfVertices(); ++nv ) {
	VertexBasePtr vertex = model()->vertex(nv);
 	AbstractFFSVertexPtr svert = 
 	  dynamic_ptr_cast<AbstractFFSVertexPtr>(vertex);
 	if(!svert) continue;
	if(vertex->getNpoint() != 3) continue;
	// check q qb allowed
	if(!vertex->isOutgoing(qk)||
	   !vertex->isOutgoing(qb)) continue;
	// check outgoing higgs allowed
	if(!vertex->isOutgoing(*ih)) continue;
  	// create the MatrixElement object
  	string objectname ("/Herwig/MatrixElements/");
  	string classname("Herwig::GeneralQQHiggs");
  	objectname += "MEPP2";
	if(iq==5) objectname += "bbbar";
	else      objectname += "ttbar";
 	objectname += (**ih).PDGName();
	GeneralQQHiggsPtr matrixElement = dynamic_ptr_cast<GeneralQQHiggsPtr>
	  (generator()->preinitCreate(classname, objectname));
	if( !matrixElement )
	  throw Exception()
	    << "QQHiggsProcessConstructor::constructDiagrams() "
	    << " Failed to construct matrix element for "
	    << qk->PDGName() << " + " << qb->PDGName() << " + " 
	    << (**ih).PDGName() << " production"
	    << Exception::runerror;
	// set the information
	matrixElement->setProcessInfo( iq, *ih, svert,_shapeOpt,
				       _process );
	// insert it
	generator()->preinitInterface(subProcess(), "MatrixElements", 
				      subProcess()->MEs().size(),
				      "insert", matrixElement->fullName());
      }
    }
  }
}
