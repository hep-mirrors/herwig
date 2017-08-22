// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HardProcessConstructor class.
//

#include "HardProcessConstructor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void HardProcessConstructor::persistentOutput(PersistentOStream & os) const {
  os << debug_ << subProcess_ << model_;
}

void HardProcessConstructor::persistentInput(PersistentIStream & is, int) {
  is >> debug_ >> subProcess_ >> model_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<HardProcessConstructor,Interfaced>
describeHerwigHardProcessConstructor("Herwig::HardProcessConstructor", "Herwig.so");

void HardProcessConstructor::Init() {

  static ClassDocumentation<HardProcessConstructor> documentation
    ("Base class for implementation of the automatic generation of hard processes");

  static Switch<HardProcessConstructor,bool> interfaceDebugME
    ("DebugME",
     "Print comparison with analytical ME",
     &HardProcessConstructor::debug_, false, false, false);
  static SwitchOption interfaceDebugMEYes
    (interfaceDebugME,
     "Yes",
     "Print the debug information",
     true);
  static SwitchOption interfaceDebugMENo
    (interfaceDebugME,
     "No",
     "Do not print the debug information",
     false);

}

void HardProcessConstructor::doinit() {
  Interfaced::doinit();
  EGPtr eg = generator();
  model_ = dynamic_ptr_cast<HwSMPtr>(eg->standardModel());
  if(!model_)
    throw InitException() << "HardProcessConstructor:: doinit() - "
			  << "The model pointer is null!"
			  << Exception::abortnow;
  if(!eg->eventHandler()) {
    throw
      InitException() << "HardProcessConstructor:: doinit() - "
		      << "The eventHandler pointer was null therefore "
		      << "could not get SubProcessHandler pointer " 
		      << Exception::abortnow;
  }
  string subProcessName = 
    eg->preinitInterface(eg->eventHandler(), "SubProcessHandlers", "get","");
  subProcess_ = eg->getObject<SubProcessHandler>(subProcessName);
  if(!subProcess_) {
    ostringstream s;
    s << "HardProcessConstructor:: doinit() - "
      << "There was an error getting the SubProcessHandler "
      << "from the current event handler. ";
    generator()->logWarning( Exception(s.str(), Exception::warning) );
  }
}

GeneralHardME::ColourStructure HardProcessConstructor::
colourFlow(const tcPDVector & extpart) const {
  PDT::Colour ina = extpart[0]->iColour();
  PDT::Colour inb = extpart[1]->iColour();
  PDT::Colour outa = extpart[2]->iColour();
  PDT::Colour outb = extpart[3]->iColour();

  // incoming colour neutral
  if(ina == PDT::Colour0 && inb == PDT::Colour0) {
    if( outa == PDT::Colour0 && outb == PDT::Colour0 ) {
      return GeneralHardME::Colour11to11;
    }
    else if( outa == PDT::Colour3 && outb == PDT::Colour3bar ) {
      return GeneralHardME::Colour11to33bar;
    } 
    else if( outa == PDT::Colour8 && outb == PDT::Colour8 ) {
      return GeneralHardME::Colour11to88;
    } 
    else
      assert(false);
  }
  // incoming 3 3
  else if(ina == PDT::Colour3 && inb == PDT::Colour3 ) {
    if( outa == PDT::Colour3 && outb == PDT::Colour3 ) {
      return GeneralHardME::Colour33to33;
    }
    else if( outa == PDT::Colour6 && outb == PDT::Colour0 ) {
      return GeneralHardME::Colour33to61;
    }
    else if( outa == PDT::Colour0 && outb == PDT::Colour6 ) {
      return GeneralHardME::Colour33to16;
    }
    else if ( outa == PDT::Colour0 && outb == PDT::Colour3bar) {
      return GeneralHardME::Colour33to13bar;
    }
    else if ( outb == PDT::Colour0 && outa == PDT::Colour3bar) {
      return GeneralHardME::Colour33to3bar1;
    }
    else if ( outa == PDT::Colour8 && outb == PDT::Colour3bar) {
      return GeneralHardME::Colour33to83bar;
    }
    else if ( outb == PDT::Colour8 && outa == PDT::Colour3bar) {
      return GeneralHardME::Colour33to3bar8;
    }
    else
      assert(false);
  }
  // incoming 3bar 3bar
  else if(ina == PDT::Colour3bar && inb == PDT::Colour3bar ) {
    if( outa == PDT::Colour3bar && outb == PDT::Colour3bar ) {
      return GeneralHardME::Colour3bar3barto3bar3bar;
    }
    else if( outa == PDT::Colour6bar && outb == PDT::Colour0) {
      return GeneralHardME::Colour3bar3barto6bar1;
    }
    else if ( outa == PDT::Colour0 && outb == PDT::Colour6bar ) {
      return GeneralHardME::Colour3bar3barto16bar;
    }
    else if ( outa == PDT::Colour0 && outb == PDT::Colour3) {
      return GeneralHardME::Colour3bar3barto13;
    }
    else if ( outb == PDT::Colour0 && outa == PDT::Colour3) {
      return GeneralHardME::Colour3bar3barto31;
    }
    else if ( outa == PDT::Colour8 && outb == PDT::Colour3) {
      return GeneralHardME::Colour3bar3barto83;
    }
    else if ( outb == PDT::Colour8 && outa == PDT::Colour3) {
      return GeneralHardME::Colour3bar3barto38;
    }
    else
      assert(false);
  }
  // incoming 3 3bar
  else if(ina == PDT::Colour3 && inb == PDT::Colour3bar ) {
    if( outa == PDT::Colour0 && outb == PDT::Colour0 ) {
      return GeneralHardME::Colour33barto11;
    }
    else if( outa == PDT::Colour3 && outb == PDT::Colour3bar ) {
      return GeneralHardME::Colour33barto33bar;
    }
    else if( outa == PDT::Colour8 && outb == PDT::Colour8 ) {
      return GeneralHardME::Colour33barto88;
    }
    else if( outa == PDT::Colour8 && outb == PDT::Colour0 ) {
      return GeneralHardME::Colour33barto81;
    }
    else if( outa == PDT::Colour0 && outb == PDT::Colour8 ) {
      return GeneralHardME::Colour33barto18;
    }
    else if( outa == PDT::Colour6 && outb == PDT::Colour6bar) {
      return GeneralHardME::Colour33barto66bar;
    }
    else if( outa == PDT::Colour6bar && outb == PDT::Colour6) {
      return GeneralHardME::Colour33barto6bar6;
    }
    else
      assert(false);
  }
  // incoming 88
  else if(ina == PDT::Colour8 && inb == PDT::Colour8 ) {
    if( outa == PDT::Colour0 && outb == PDT::Colour0 ) {
      return GeneralHardME::Colour88to11;
    }
    else if( outa == PDT::Colour3 && outb == PDT::Colour3bar ) {
      return GeneralHardME::Colour88to33bar;
    }
    else if( outa == PDT::Colour8 && outb == PDT::Colour8 ) {
      return GeneralHardME::Colour88to88;
    }
    else if( outa == PDT::Colour8 && outb == PDT::Colour0 ) {
      return GeneralHardME::Colour88to81;
    }
    else if( outa == PDT::Colour0 && outb == PDT::Colour8 ) {
      return GeneralHardME::Colour88to18;
    }
    else if( outa == PDT::Colour6 && outb == PDT::Colour6bar ) {
      return GeneralHardME::Colour88to66bar;
    }    
    else
      assert(false);
  }
  // incoming 38
  else if(ina == PDT::Colour3 && inb == PDT::Colour8 ) {
    if(outa == PDT::Colour3 && outb == PDT::Colour0) {
      return GeneralHardME::Colour38to31;
    }
    else if(outa == PDT::Colour0 && outb == PDT::Colour3) {
      return GeneralHardME::Colour38to13;
    }
    else if(outa == PDT::Colour3 && outb == PDT::Colour8) {
      return GeneralHardME::Colour38to38;
    }
    else if(outa == PDT::Colour8 && outb == PDT::Colour3) {
      return GeneralHardME::Colour38to83;
    }
    else if(outa == PDT::Colour3bar && outb == PDT::Colour6){
      return GeneralHardME::Colour38to3bar6;
    }
    else if(outa == PDT::Colour6 && outb == PDT::Colour3bar) {
      return GeneralHardME::Colour38to63bar;
    }
    else if(outa == PDT::Colour3bar && outb == PDT::Colour3bar) {
      return GeneralHardME::Colour38to3bar3bar;
    }
    else
      assert(false);
  }
  // incoming 3bar8
  else if(ina == PDT::Colour3bar && inb == PDT::Colour8 ) {   
    if(outa == PDT::Colour3bar && outb == PDT::Colour0 ) {
      return GeneralHardME::Colour3bar8to3bar1;
    }
    else if(outa == PDT::Colour0 && outb == PDT::Colour3bar) {
      return GeneralHardME::Colour3bar8to13bar;
    }
    else if(outa == PDT::Colour3bar && outb == PDT::Colour8 ) {
      return GeneralHardME::Colour3bar8to3bar8;
    }
    else if(outa == PDT::Colour8 && outb == PDT::Colour3bar) {
      return GeneralHardME::Colour3bar8to83bar;
    }
    else if(outa == PDT::Colour3 && outb == PDT::Colour3) {
      return GeneralHardME::Colour3bar8to33;
    }
    else
      assert(false);
  }
  // unknown colour flow
  else 
    assert(false);
  return GeneralHardME::UNDEFINED;
}


void HardProcessConstructor::fixFSOrder(HPDiagram & diag) {
  tcPDPtr psa = getParticleData(diag.incoming.first);
  tcPDPtr psb = getParticleData(diag.incoming.second);
  tcPDPtr psc = getParticleData(diag.outgoing.first);
  tcPDPtr psd = getParticleData(diag.outgoing.second);

  //fix a spin order
  if( psc->iSpin() < psd->iSpin() ) {
    swap(diag.outgoing.first, diag.outgoing.second);
    if(diag.channelType == HPDiagram::tChannel) {
      diag.ordered.second = !diag.ordered.second;
    }
    return;
  }
  
  if( psc->iSpin() == psd->iSpin() && 
      psc->id() < 0 && psd->id() > 0 ) {
    swap(diag.outgoing.first, diag.outgoing.second);
    if(diag.channelType == HPDiagram::tChannel) {
      diag.ordered.second = !diag.ordered.second;
    }
    return;
  }
}

void HardProcessConstructor::assignToCF(HPDiagram & diag) {
  if(diag.channelType == HPDiagram::tChannel) {
    if(diag.ordered.second) tChannelCF(diag);
    else                    uChannelCF(diag);
  }
  else if(diag.channelType == HPDiagram::sChannel) {
    sChannelCF(diag);
  }
  else if (diag.channelType == HPDiagram::fourPoint) {
    fourPointCF(diag);
  }
  else 
    assert(false);
}

void HardProcessConstructor::tChannelCF(HPDiagram & diag) {
  tcPDPtr ia = getParticleData(diag.incoming.first );
  tcPDPtr ib = getParticleData(diag.incoming.second);
  tcPDPtr oa = getParticleData(diag.outgoing.first );
  tcPDPtr ob = getParticleData(diag.outgoing.second);
  PDT::Colour ina  = ia->iColour();
  PDT::Colour inb  = ib->iColour();
  PDT::Colour outa = oa->iColour();
  PDT::Colour outb = ob->iColour();
  vector<CFPair> cfv(1, make_pair(0, 1.));
  if(diag.intermediate->iColour() == PDT::Colour0) {
    if(ina==PDT::Colour0) {
      cfv[0] = make_pair(0, 1);
    }
    else if(ina==PDT::Colour3 || ina==PDT::Colour3bar) {
      if( inb == PDT::Colour0 ) {
	cfv[0] = make_pair(0, 1);
      }
      else if(inb==PDT::Colour3 || outb==PDT::Colour3bar) {
	cfv[0] = make_pair(2, 1);
      }
      else if(inb==PDT::Colour8) {
	cfv[0] = make_pair(2, 1);
      }
    }
    else if(ina==PDT::Colour8) {
      if( inb == PDT::Colour0 ) {
	cfv[0] = make_pair(0, 1);
      }
      else if(inb==PDT::Colour3 || outb==PDT::Colour3bar) {
	cfv[0] = make_pair(2, 1);
      }
      else if(inb==PDT::Colour8) {
	cfv[0] = make_pair(7, 1);
      }
    }
  }
  else if(diag.intermediate->iColour() == PDT::Colour8) {
    if(ina==PDT::Colour8&&outa==PDT::Colour8&&
       inb==PDT::Colour8&&outb==PDT::Colour8) {
      cfv[0]=make_pair(2,  2.);
      cfv.push_back(make_pair(3, -2.));
      cfv.push_back(make_pair(1, -2.));
      cfv.push_back(make_pair(4,  2.));
    }
    else if(ina==PDT::Colour8&&outa==PDT::Colour0&&
	    inb==PDT::Colour8&&outb==PDT::Colour8&&
	    oa->iSpin()==PDT::Spin0) {
      cfv[0] = make_pair(0,-1);
    }
    else if(ina==PDT::Colour8&&outa==PDT::Colour8&&
	    inb==PDT::Colour8&&outb==PDT::Colour0&&
	    ob->iSpin()==PDT::Spin0) {
      cfv[0] = make_pair(0,-1);
    }
  } 
  else if(diag.intermediate->iColour() == PDT::Colour3 ||
	  diag.intermediate->iColour() == PDT::Colour3bar) {
    if(outa == PDT::Colour0 || outb == PDT::Colour0) {
      if( outa != PDT::Colour6    && outb != PDT::Colour6   &&
	  outa != PDT::Colour6bar && outb != PDT::Colour6bar) {
	cfv[0] = make_pair(0,1.);
      }
      else {
	cfv[0] = make_pair(0,0.5);
	cfv.push_back(make_pair(1,0.5));
      }
    }
    else if(outa==PDT::Colour6 && outb==PDT::Colour3bar) {
      cfv[0] = make_pair(4,1.);
      cfv.push_back(make_pair(5,1.));
    }
    else if(outa==PDT::Colour6 && outb==PDT::Colour6bar) {
      cfv[0] = make_pair(4, 1.);
      for(unsigned int ix=5;ix<8;++ix)
	cfv.push_back(make_pair(ix,1.));
    }
    else if(outa==PDT::Colour6 || outa ==PDT::Colour6bar ||
	    outb==PDT::Colour6 || outb ==PDT::Colour6bar ) {
      assert(false);
    }
    else if(ina==PDT::Colour3    && inb==PDT::Colour3    ) {
      if((outa==PDT::Colour0 && outb==PDT::Colour3bar)||
	 (outb==PDT::Colour0 && outa==PDT::Colour3bar))
	cfv[0] = make_pair(0,1.);
      else if((outa==PDT::Colour8 && outb==PDT::Colour3bar)||
	      (outb==PDT::Colour8 && outa==PDT::Colour3bar))
	cfv[0] = make_pair(1,1.);
    }
    else if(ina==PDT::Colour3bar && inb==PDT::Colour3bar ) {
      if((outa==PDT::Colour0 && outb==PDT::Colour3)||
	 (outb==PDT::Colour0 && outa==PDT::Colour3))
	cfv[0] = make_pair(0,1.);
      else if((outa==PDT::Colour8 && outb==PDT::Colour3)||
	      (outb==PDT::Colour8 && outa==PDT::Colour3))
	cfv[0] = make_pair(1,1.);
    }
    else if((ina==PDT::Colour3    && inb==PDT::Colour8) ||
	    (ina==PDT::Colour3bar && inb==PDT::Colour8) ||
	    (inb==PDT::Colour3    && ina==PDT::Colour8) ||
	    (inb==PDT::Colour3bar && ina==PDT::Colour8) ) {
      if((outa==PDT::Colour3    && outb==PDT::Colour3    ) ||
	 (outa==PDT::Colour3bar && outb==PDT::Colour3bar)) {
	cfv[0] = make_pair(1,1.);
      }
    }
  }
  else if(diag.intermediate->iColour() == PDT::Colour6 ||
	  diag.intermediate->iColour() == PDT::Colour6bar) {
    if(ina==PDT::Colour8 && inb==PDT::Colour8) {
      cfv[0] = make_pair(0, 1.);
      for(unsigned int ix=1;ix<4;++ix)
	cfv.push_back(make_pair(ix,1.));
      for(unsigned int ix=4;ix<8;++ix)
	cfv.push_back(make_pair(ix,1.));
    }
    else if(outa==PDT::Colour3bar && outb==PDT::Colour6) {
      cfv[0] = make_pair(0,1.);
      for(unsigned int ix=1;ix<4;++ix)
	cfv.push_back(make_pair(ix,1.));
    }
    else if(outa==PDT::Colour6 && outb==PDT::Colour3bar) {
      cfv[0] = make_pair(4,1.);
      cfv.push_back(make_pair(5,1.));
    }
  }
  diag.colourFlow = cfv;
}
 
void HardProcessConstructor::uChannelCF(HPDiagram & diag) {
  PDT::Colour offshell = diag.intermediate->iColour();
  PDT::Colour ina  = getParticleData(diag.incoming.first )->iColour();
  PDT::Colour inb  = getParticleData(diag.incoming.second)->iColour();
  PDT::Colour outa = getParticleData(diag.outgoing.first )->iColour();
  PDT::Colour outb = getParticleData(diag.outgoing.second)->iColour();
  vector<CFPair> cfv(1, make_pair(1, 1.));
  if(offshell == PDT::Colour8) {
    if(outa == PDT::Colour0 &&
       outb == PDT::Colour0) {
      cfv[0].first = 0;
    }
    else if( outa != outb ) {
      if(outa == PDT::Colour0 || 
	 outb == PDT::Colour0) {
	cfv[0].first = 0;
      }
      else if(ina  == PDT::Colour3 && inb  == PDT::Colour8 &&
	      outb == PDT::Colour3 && outa == PDT::Colour8) {
	tPDPtr off = diag.intermediate;
	if(off->CC()) off=off->CC();
	if(off->iSpin()!=PDT::Spin1Half ||
	   diag.vertices.second->allowed(off->id(),diag.outgoing.first,diag.incoming.second)) {
	  cfv[0].first = 0;
	  cfv.push_back(make_pair(1, -1.));
	}
	else {
	  cfv[0].first = 1;
	  cfv.push_back(make_pair(0, -1.));
	}
      }
      else if(ina  == PDT::Colour3bar && inb  == PDT::Colour8 &&
	      outb == PDT::Colour3bar && outa == PDT::Colour8) {
	tPDPtr off = diag.intermediate;
	if(off->CC()) off=off->CC();
	if(off->iSpin()!=PDT::Spin1Half ||
	   diag.vertices.second->allowed(diag.outgoing.first,off->id(),diag.incoming.second)) {
	  cfv[0].first = 0;
	  cfv.push_back(make_pair(1, -1.));
	}
	else {
	  cfv[0].first = 1;
	  cfv.push_back(make_pair(0, -1.));
	}
      }
      else {
	cfv[0].first = 0;
	cfv.push_back(make_pair(1, -1.));
      }
    }
    else if(outa==PDT::Colour8&&ina==PDT::Colour8) {
      cfv[0]=make_pair(4, 2.);
      cfv.push_back(make_pair(5, -2.));
      cfv.push_back(make_pair(0, -2.));
      cfv.push_back(make_pair(2,  2.));
    }
  }
  else if(offshell == PDT::Colour3 || offshell == PDT::Colour3bar) {
    if( outa == PDT::Colour0 || outb == PDT::Colour0 ) {
      if( outa != PDT::Colour6    && outb != PDT::Colour6   &&
	  outa != PDT::Colour6bar && outb != PDT::Colour6bar) {
	cfv[0] = make_pair(0,1.);
      }
      else {
	cfv[0] = make_pair(0,0.5);
	cfv.push_back(make_pair(1,0.5));
      }
    }
    else if(outa==PDT::Colour3bar && outb==PDT::Colour6) {
      cfv[0] = make_pair(4,1.);
      cfv.push_back(make_pair(5,1.));
    }
    else if(outa==PDT::Colour6 && outb==PDT::Colour3bar) {
      cfv[0] = make_pair(0,1.);
      for(int ix=0; ix<4;++ix)
	cfv.push_back(make_pair(ix,1.));
    }
    else if(outa==PDT::Colour6bar && outb==PDT::Colour6) {
      cfv[0] = make_pair(4,1.);
      for(int ix=5; ix<8;++ix)
	cfv.push_back(make_pair(ix,1.));
    }
    else if(ina==PDT::Colour0 && inb==PDT::Colour0) {
      cfv[0] = make_pair(0,1.);
    }
    else if(ina==PDT::Colour3    && inb==PDT::Colour3    ) {
      if((outa==PDT::Colour0 && outb==PDT::Colour3bar)||
	 (outb==PDT::Colour0 && outa==PDT::Colour3bar))
	cfv[0] = make_pair(0,1.);
      else if((outa==PDT::Colour8 && outb==PDT::Colour3bar)||
	      (outb==PDT::Colour8 && outa==PDT::Colour3bar))
	cfv[0] = make_pair(2,1.);
    }
    else if(ina==PDT::Colour3bar && inb==PDT::Colour3bar ) {
      if((outa==PDT::Colour0 && outb==PDT::Colour3)||
	 (outb==PDT::Colour0 && outa==PDT::Colour3))
	cfv[0] = make_pair(0,1.);
      else if((outa==PDT::Colour8 && outb==PDT::Colour3)||
	      (outb==PDT::Colour8 && outa==PDT::Colour3))
	cfv[0] = make_pair(2,1.);
    }
    else if(((ina==PDT::Colour3    && inb==PDT::Colour8) ||
	     (ina==PDT::Colour3bar && inb==PDT::Colour8) ||
	     (inb==PDT::Colour3    && ina==PDT::Colour8) ||
	     (inb==PDT::Colour3bar && ina==PDT::Colour8)) &&
	    ((outa==PDT::Colour3    && outb==PDT::Colour3    ) ||
	     (outa==PDT::Colour3bar && outb==PDT::Colour3bar))) {
      cfv[0] = make_pair(2, 1.);
    }
    else if(( ina==PDT::Colour3    &&  inb==PDT::Colour3bar && 
	      outa==PDT::Colour3    && outb==PDT::Colour3bar)) {
      cfv[0] = make_pair(2, 1.);
      cfv.push_back(make_pair(3,-1.));
    }
  }
  else if( offshell == PDT::Colour0 ) {
    if(ina==PDT::Colour0) {
      cfv[0] = make_pair(0, 1);
    }
    else if(ina==PDT::Colour3 || ina==PDT::Colour3bar) {
      if( inb == PDT::Colour0 ) {
	cfv[0] = make_pair(0, 1);
      }
      else if(inb==PDT::Colour3 || outb==PDT::Colour3bar) {
	cfv[0] = make_pair(3, 1);
      }
      else if(inb==PDT::Colour8) {
	cfv[0] = make_pair(2, 1);
      }
    }
    else if(ina==PDT::Colour8) {
      if( inb == PDT::Colour0 ) {
	cfv[0] = make_pair(0, 1);
      }
      else if(inb==PDT::Colour3 || outb==PDT::Colour3bar) {
	cfv[0] = make_pair(2, 1);
      }
      else if(inb==PDT::Colour8) {
	cfv[0] = make_pair(8, 1);
      }
    }
  }
  else if(diag.intermediate->iColour() == PDT::Colour6 ||
	  diag.intermediate->iColour() == PDT::Colour6bar) {
    if(ina==PDT::Colour8 && inb==PDT::Colour8) {
      cfv[0] = make_pair(0, 1.);
      for(unsigned int ix=1;ix<4;++ix)
	cfv.push_back(make_pair(ix,1.));
      for(unsigned int ix=8;ix<12;++ix)
	cfv.push_back(make_pair(ix,1.));
    }
    else if(outa==PDT::Colour3bar && outa==PDT::Colour6) {
      cfv[0] = make_pair(4, 1.);
      cfv.push_back(make_pair(5,1.));
    }
    else if(outa==PDT::Colour6 && outa==PDT::Colour3bar) {
      cfv[0] = make_pair(0, 1.);
      for(unsigned int ix=1;ix<4;++ix)
	cfv.push_back(make_pair(ix,1.));
    }
  }
  diag.colourFlow = cfv;
}

void HardProcessConstructor::sChannelCF(HPDiagram & diag) {
  tcPDPtr pa = getParticleData(diag.incoming.first);
  tcPDPtr pb = getParticleData(diag.incoming.second);
  PDT::Colour ina = pa->iColour();
  PDT::Colour inb = pb->iColour();
  PDT::Colour offshell = diag.intermediate->iColour();
  tcPDPtr pc = getParticleData(diag.outgoing.first);
  tcPDPtr pd = getParticleData(diag.outgoing.second);
  PDT::Colour outa = pc->iColour();
  PDT::Colour outb = pd->iColour();
  vector<CFPair> cfv(1);
  if(offshell == PDT::Colour8) {
    if(ina  == PDT::Colour0 || inb  == PDT::Colour0 || 
       outa == PDT::Colour0 || outb == PDT::Colour0) {
      cfv[0] = make_pair(0, 1);
    }
    else {
      bool incol   = ina  == PDT::Colour8 && inb  == PDT::Colour8;
      bool outcol  = outa == PDT::Colour8 && outb == PDT::Colour8;
      bool intrip  = ina  == PDT::Colour3 && inb  == PDT::Colour3bar;
      bool outtrip = outa == PDT::Colour3 && outb == PDT::Colour3bar;
      bool outsex  = outa == PDT::Colour6 && outb == PDT::Colour6bar;
      bool outsexb = outa == PDT::Colour6bar && outb == PDT::Colour6;
      if(incol || outcol) {
	// Require an additional minus sign for a scalar/fermion
	// 33bar final state due to the way the vertex rules are defined.
	int prefact(1);
	if( ((pc->iSpin() == PDT::Spin1Half && pd->iSpin() == PDT::Spin1Half) ||
	     (pc->iSpin() == PDT::Spin0     && pd->iSpin() == PDT::Spin0    )) &&
	    (outa        == PDT::Colour3   && outb        == PDT::Colour3bar) )
	  prefact = -1;
	if(incol && outcol) {
	  cfv[0] = make_pair(0, -2.);
	  cfv.push_back(make_pair(1,  2.));
	  cfv.push_back(make_pair(3,  2.));
	  cfv.push_back(make_pair(5,  -2.));
	}
	else if(incol && outsex) {
	  cfv[0].first = 4;
	  cfv[0].second =  prefact;
	  for(unsigned int ix=1;ix<4;++ix)
	    cfv.push_back(make_pair(4+ix, prefact));
	  for(unsigned int ix=0;ix<4;++ix)
	    cfv.push_back(make_pair(8+ix,-prefact));
	}
	else {
	  cfv[0].first = 0;
	  cfv[0].second = -prefact;
	  cfv.push_back(make_pair(1, prefact));
	}
      }
      else if( (  intrip && !outtrip ) || 
	       ( !intrip &&  outtrip ) ) {
	if(!outsex)
	  cfv[0] = make_pair(0, 1);
	else {
	  cfv[0] = make_pair(0, 1.);
	  for(unsigned int ix=0;ix<3;++ix)
	    cfv.push_back(make_pair(ix+1, 1.));
	}
      }
      else if((intrip && outsex) || (intrip && outsexb)) {
	cfv[0] = make_pair(0,1.);
	for(int ix=1; ix<4; ++ix)
	  cfv.push_back(make_pair(ix,1.));
      }
      else
	cfv[0] = make_pair(1, 1);
    }
  }
  else if(offshell == PDT::Colour0) {
    if( ina == PDT::Colour0 ) {
      cfv[0] = make_pair(0, 1);
    }
    else if(ina==PDT::Colour3 || ina==PDT::Colour3bar) {
      if( outa == PDT::Colour0 ) {
	cfv[0] = make_pair(0, 1);
      }
      else if(outa==PDT::Colour3 || outa==PDT::Colour3bar) {
	cfv[0] = make_pair(3, 1);
      }
      else if(outa==PDT::Colour8) {
	cfv[0] = make_pair(2, 1);
      }
      else if(outa==PDT::Colour6 || outa==PDT::Colour6bar) {
	cfv[0] = make_pair(8, 1.);
	cfv.push_back(make_pair(9,1.));
      }
      else
	assert(false);
    }
    else if(ina==PDT::Colour8) {
      if( outa == PDT::Colour0 ) {
	cfv[0] = make_pair(0, 1);
      }
      else if(outa==PDT::Colour3 || outb==PDT::Colour3bar) {
	cfv[0] = make_pair(2, 1);
      }
      else if(outa==PDT::Colour8) {
	cfv[0] = make_pair(6, 1);
      }
    }
  }
  else if(offshell == PDT::Colour3 || offshell == PDT::Colour3bar) {
    if(outa == PDT::Colour6    || outa == PDT::Colour6bar || 
       outb == PDT::Colour6bar || outb == PDT::Colour6) {
      cfv[0] = make_pair(6, 1.);
      cfv.push_back(make_pair(7,1.));
    }
    else if((ina  == PDT::Colour3    && inb  == PDT::Colour3) ||
	    (ina  == PDT::Colour3bar && inb  == PDT::Colour3bar)) {
      if((outa == PDT::Colour3    && outb == PDT::Colour3   ) ||
	 (outa == PDT::Colour3bar && outb == PDT::Colour3bar)) {
	cfv[0]      = make_pair(2, 1.);
	cfv.push_back(make_pair(3,-1.));
      }
      else
	cfv[0] = make_pair(0,1.);
    }
    else if(((ina==PDT::Colour3    && inb==PDT::Colour8) ||
	     (ina==PDT::Colour3bar && inb==PDT::Colour8) ||
	     (inb==PDT::Colour3    && ina==PDT::Colour8) ||
	     (inb==PDT::Colour3bar && ina==PDT::Colour8) ) &&
	    ((outa==PDT::Colour3    && outb==PDT::Colour3    ) ||
	     (outa==PDT::Colour3bar && outb==PDT::Colour3bar))) {
      cfv[0] = make_pair(0,1.);
    }
    else {
      if(outa == PDT::Colour0 || outb == PDT::Colour0)
	cfv[0] = make_pair(0, 1);
      else
	cfv[0] = make_pair(1, 1);
    }
  }
  else if( offshell == PDT::Colour6 || offshell == PDT::Colour6bar) {
    if((ina  == PDT::Colour3    && inb  == PDT::Colour3    &&
	outa == PDT::Colour3    && outb == PDT::Colour3   ) ||
       (ina  == PDT::Colour3bar && inb  == PDT::Colour3bar &&
	outa == PDT::Colour3bar && outb == PDT::Colour3bar)) {
      cfv[0]      = make_pair(2,0.5);
      cfv.push_back(make_pair(3,0.5));
    }
    else if((ina  == PDT::Colour3    && inb  == PDT::Colour3    &&
	     ((outa == PDT::Colour6    && outb == PDT::Colour0)||
	      (outb == PDT::Colour6    && outa == PDT::Colour0))) ||
	    (ina  == PDT::Colour3bar && inb  == PDT::Colour3bar &&
	     ((outa == PDT::Colour6bar && outb == PDT::Colour0)||
	      (outb == PDT::Colour6bar && outa == PDT::Colour0)))) {
      cfv[0]      = make_pair(0,0.5);
      cfv.push_back(make_pair(1,0.5));
    }
    else
      assert(false);
  }
  else {
    if(outa == PDT::Colour0 || outb == PDT::Colour0)
      cfv[0] = make_pair(0, 1);
    else
      cfv[0] = make_pair(1, 1);
  }  
  diag.colourFlow = cfv; 
}

void HardProcessConstructor::fourPointCF(HPDiagram & diag) {
  // count the colours
  unsigned int noct(0),ntri(0),nsng(0),nsex(0);
  for(unsigned int ix=0;ix<4;++ix) {
    PDT::Colour col = getParticleData(diag.ids[ix])->iColour();
    if(col==PDT::Colour0)                            ++nsng;
    else if(col==PDT::Colour3||col==PDT::Colour3bar) ++ntri;
    else if(col==PDT::Colour8)                       ++noct;
    else if(col==PDT::Colour6||col==PDT::Colour6bar) ++nsex;
  }
  if(nsng==4 || (ntri==2&&nsng==2) || 
     (noct==3            && nsng==1) ||
     (ntri==2 && noct==1 && nsng==1) ) {
    vector<CFPair> cfv(1,make_pair(0,1));
    diag.colourFlow = cfv;
  }
  else if(noct==4) {
    // flows for SSVV, VVVV is handled in me class
    vector<CFPair> cfv(6);
    cfv[0] = make_pair(0, -2.);
    cfv[1] = make_pair(1, -2.);
    cfv[2] = make_pair(2, +4.);
    cfv[3] = make_pair(3, -2.);
    cfv[4] = make_pair(4, +4.);
    cfv[5] = make_pair(5, -2.);
    diag.colourFlow = cfv;
  }
  else if(ntri==2&&noct==2) {
    vector<CFPair> cfv(2);
    cfv[0] = make_pair(0, 1);
    cfv[1] = make_pair(1, 1);
    diag.colourFlow = cfv;
  }
  else if(nsex==2&&noct==2) {
    vector<CFPair> cfv;
    for(unsigned int ix=0;ix<4;++ix)
      cfv.push_back(make_pair(ix  ,2.));
    for(unsigned int ix=0;ix<8;++ix)
      cfv.push_back(make_pair(4+ix,1.));
    diag.colourFlow = cfv;
  }
  else
    assert(false);
}

namespace {
  // Helper functor for find_if in duplicate function.
  class SameDiagramAs {
  public:
    SameDiagramAs(const HPDiagram & diag) : a(diag) {}
    bool operator()(const HPDiagram & b) const {
      return a == b;
    }
  private:
    HPDiagram a;
  };
}

bool HardProcessConstructor::duplicate(const HPDiagram & diag, 
				       const HPDVector & group) const {
  //find if a duplicate diagram exists
  HPDVector::const_iterator it = 
    find_if(group.begin(), group.end(), SameDiagramAs(diag));
  return it != group.end();
} 
