// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralFourBodyDecayer class.
//

#include "GeneralFourBodyDecayer.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void GeneralFourBodyDecayer::persistentOutput(PersistentOStream & os) const {
  os << _incoming << _outgoing << _diagrams << _reftag << _reftagcc
     << _widthopt << _colour << _colourLargeNC << _nflow << _diagmap;
}

void GeneralFourBodyDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incoming >> _outgoing >> _diagrams >> _reftag >> _reftagcc
     >> _widthopt >> _colour >> _colourLargeNC >> _nflow >> _diagmap;
}

DescribeAbstractClass<GeneralFourBodyDecayer,DecayIntegrator>
describeGeneralFourBodyDecayer("Herwig::GeneralFourBodyDecayer",
			       "Herwig.so");

void GeneralFourBodyDecayer::Init() {

  static ClassDocumentation<GeneralFourBodyDecayer> documentation
    ("The GeneralFourBodyDecayer class is the base class for the implementation "
     "of all four-body decays based on spin structures in Herwig++.");

  static Switch<GeneralFourBodyDecayer,unsigned int> interfaceWidthOption
    ("WidthOption",
     "Option for the treatment of the widths of the intermediates",
     &GeneralFourBodyDecayer::_widthopt, 1, false, false);
  static SwitchOption interfaceWidthOptionFixed
    (interfaceWidthOption,
     "Fixed",
     "Use fixed widths",
     1);
  static SwitchOption interfaceWidthOptionRunning
    (interfaceWidthOption,
     "Running",
     "Use running widths",
     2);
  static SwitchOption interfaceWidthOptionZero
    (interfaceWidthOption,
     "Zero",
     "Set the widths to zero",
     3);

}

ParticleVector GeneralFourBodyDecayer::decay(const Particle & parent,
					     const tPDVector & children) const {
  // return empty vector if products heavier than parent
  Energy mout(ZERO);
  for(tPDVector::const_iterator it=children.begin();
      it!=children.end();++it) mout+=(**it).massMin();
  if(mout>parent.mass()) return ParticleVector();
  // generate the decay
  bool cc;
  int imode=modeNumber(cc,parent.dataPtr(),children);
  // generate the kinematics
  ParticleVector decay=generate(generateIntermediates(),cc,imode,parent);
  // make the colour connections
  colourConnections(parent, decay);
  // return the answer
  return decay;
}

int GeneralFourBodyDecayer::modeNumber(bool & cc, tcPDPtr parent,
				       const tPDVector & children) const {
  assert( !_reftag.empty() && !_reftagcc.empty() );
  // check number of outgoing particles
  if( children.size() != 4 || abs(parent->id()) != _incoming->id() ) return -1;
  OrderedParticles testmode(children.begin(), children.end());
  OrderedParticles::const_iterator dit = testmode.begin();
  string testtag(parent->name() + "->");
  for( unsigned int i = 1; dit != testmode.end(); ++dit, ++i) {
    testtag += (**dit).name();
    if( i != 3 ) testtag += string(",");
  }
  if( testtag == _reftag ) {
    cc = false;
    return 0;
  }
  else if ( testtag == _reftagcc ) {
    cc = true;
    return 0;
  }
  else return -1;
}
 
bool GeneralFourBodyDecayer::setDecayInfo(PDPtr incoming,
					  vector<PDPtr> outgoing,
					  const vector<NBDiagram> & process) {
  // set the member variables from the info supplied
  // external particles
  _incoming        = incoming;
  _outgoing        = outgoing;
  _diagrams        = process;
  assert( _outgoing.size() == 4 );
  // Construct reference tags for testing in modeNumber function
  OrderedParticles refmode(_outgoing.begin(), _outgoing.end());
  OrderedParticles::const_iterator dit = refmode.begin();
  _reftag = _incoming->name() + "->";
  for( ; dit != refmode.end(); ++dit) {
    if( dit != refmode.begin() )  _reftag += string(",");
    _reftag += (**dit).name();
  }
  //CC-mode
  refmode.clear();
  _reftagcc = _incoming->CC() ? _incoming->CC()->name() : _incoming->name();
  _reftagcc += "->";
  for( unsigned int i = 0;  i < _outgoing.size(); ++i ) {
    if( _outgoing[i]->CC() ) refmode.insert( _outgoing[i]->CC() );
    else refmode.insert( _outgoing[i] );
  }
  dit = refmode.begin();
  for( ; dit != refmode.end(); ++dit) {
    if( dit != refmode.begin() )  _reftag += string(",");
    _reftagcc += (**dit).name();
  }
  // set the colour factors and return the answer
  return setColourFactors();
}

Energy GeneralFourBodyDecayer::partialWidth(tPDPtr inpart,
					    OrderedParticles outgoing) const {
  tPDVector temp = tPDVector(outgoing.begin(),outgoing.end());
  bool cc=false;
  int imode = modeNumber(cc,inpart,temp);
  if(imode<0) return ZERO;
  else return initializePhaseSpaceMode(0,true);  
}

namespace {
  bool addIntermediate(DecayPhaseSpaceChannelPtr newChannel,
		       const NBVertex &vertex,int &order,
		       int & ioff, int & loc) {
    if(vertex.vertices.size()!=2) return false;
    assert(!vertex.vertices.begin()->second.incoming);
    int first = order/pow(10,loc);
    order = order%int(pow(10,loc));
    --loc;
    // one off-shell
    if((++vertex.vertices.begin())->second.incoming) {
      newChannel->addIntermediate(vertex.incoming,0,0.0,ioff,first);
      --ioff;
      return addIntermediate(newChannel,(++vertex.vertices.begin())->second,
			     order,ioff,loc);
    }
    // no off-shell
    else {
      int second = order/pow(10,loc);
      order = order%int(pow(10,loc));
      --loc;
      newChannel->addIntermediate(vertex.incoming,0,0.0,first,second);
      --ioff;
      return true;
    }
  }
}

void GeneralFourBodyDecayer::doinit() {
  DecayIntegrator::doinit();
  // create the phase space integrator
  tPDVector extpart(1,_incoming);
  extpart.insert(extpart.end(),_outgoing.begin(),_outgoing.end());
  // create the integration channels for the decay
  DecayPhaseSpaceModePtr mode(new_ptr(DecayPhaseSpaceMode(extpart,this,true)));
  DecayPhaseSpaceChannelPtr newChannel;
  // create the phase-space channels for the integration
  unsigned int nmode(0);
  for(unsigned int ix=0;ix<_diagrams.size();++ix) {
    // check not four point
    if(_diagrams[ix].vertices.size()!=2) continue;
    // create the new channel
    newChannel=new_ptr(DecayPhaseSpaceChannel(mode));
    int order = _diagrams[ix].channelType;
    int ioff = -1;
    int loc=3;
    // two off-shell particles
    if(_diagrams[ix].vertices.begin()->second.incoming) {
      newChannel->addIntermediate(extpart[0],0,0.0,-1,-2);
      if(!addIntermediate(newChannel,(  _diagrams[ix].vertices.begin())->second,
			  order,ioff,loc)) continue;
      if(!addIntermediate(newChannel,(++_diagrams[ix].vertices.begin())->second,
			  order,ioff,loc)) continue;
    }
    // only one off-shell particle
    else {
      int first = order/pow(10,loc);
      order = order%int(pow(10,loc));
      --loc;
      newChannel->addIntermediate(extpart[0],0,0.0,-1,first);
      --ioff;
      if(!addIntermediate(newChannel,(++_diagrams[ix].vertices.begin())->second,
			  order,ioff,loc)) continue;
    }
    _diagmap.push_back(ix);
    mode->addChannel(newChannel);
    cerr << "testing channel " << *newChannel << "\n";
    ++nmode;
  }
  if(nmode==0) {
    string mode = extpart[0]->PDGName() + "->";
    for(unsigned int ix=1;ix<extpart.size();++ix) mode += extpart[ix]->PDGName() + " ";
    throw Exception() << "No decay channels in GeneralFourBodyDecayer::"
		      << "doinit() for " << mode << "\n" << Exception::runerror;
  }
  // add the mode
  vector<double> wgt(nmode,1./double(nmode));
  addMode(mode,1.,wgt);
}

void GeneralFourBodyDecayer::colourConnections(const Particle & parent,
					       const ParticleVector & out) const {
  // first extract the outgoing particles and intermediate
  ParticleVector inter;
  ParticleVector outgoing;
  if(!generateIntermediates()) {
    outgoing=out;
  }
  else {
    // find the diagram
    unsigned int idiag = diagramMap()[mode(imode())->selectedChannel()];
//     PPtr child;
//     for(unsigned int ix=0;ix<out.size();++ix) {
//       if(out[ix]->children().empty()) child = out[ix];
//       else                            inter = out[ix];
//     }
    outgoing.resize(4);
//     switch(_diagrams[idiag].channelType) {
//     case TBDiagram::channel23:
//       outgoing[0] = child;
//       outgoing[1] = inter->children()[0];
//       outgoing[2] = inter->children()[1];
//       break;
//     case TBDiagram::channel13:
//       outgoing[0] = inter->children()[0];
//       outgoing[1] = child;
//       outgoing[2] = inter->children()[1];
//       break;
//     case TBDiagram::channel12:
//       outgoing[0] = inter->children()[0];
//       outgoing[1] = inter->children()[1];
//       outgoing[2] = child;
//       break;
//     default:
//       throw Exception() << "unknown diagram type in GeneralFourBodyDecayer::"
// 			<< "colourConnections()" << Exception::runerror;
//     }
  }
//   // extract colour of the incoming and outgoing particles
//   PDT::Colour inColour(parent.data().iColour());
//   vector<PDT::Colour> outColour;
//   vector<int> singlet,octet,triplet,antitriplet;
//   for(unsigned int ix=0;ix<outgoing.size();++ix) {
//     outColour.push_back(outgoing[ix]->data().iColour());
//     switch(outColour.back()) {
//     case PDT::Colour0   :     
//       singlet.push_back(ix);
//       break;
//     case PDT::Colour3   :     
//       triplet.push_back(ix);
//       break;
//     case PDT::Colour3bar: 
//       antitriplet.push_back(ix);
//       break;
//     case PDT::Colour8   :     
//       octet.push_back(ix);
//       break;
//     default:
//       throw Exception() << "Unknown colour for particle in GeneralFourBodyDecayer::"
// 			<< "colourConnections()" << Exception::runerror;
//     }
//   }
//   // colour neutral decaying particle
//   if     ( inColour == PDT::Colour0) {
//     // options are all neutral or triplet/antitriplet+ neutral
//     if(singlet.size()==3) return;
//     else if(singlet.size()==1&&triplet.size()==1&&antitriplet.size()==1) {
//       outgoing[triplet[0]]->antiColourNeighbour(outgoing[antitriplet[0]]);
//       // add intermediate if needed
//       if(inter&&inter->coloured()) {
// 	if(inter->dataPtr()->iColour()==PDT::Colour3)
// 	  outgoing[triplet[0]]->colourLine()->addColoured(inter);
// 	else if(inter->dataPtr()->iColour()==PDT::Colour3bar)
// 	  outgoing[triplet[0]]->colourLine()->addAntiColoured(inter);
//       }
//     }
//     else if(octet.size()==1&&triplet.size()==1&&antitriplet.size()==1) {
//       outgoing[    triplet[0]]->antiColourNeighbour(outgoing[octet[0]]);
//       outgoing[antitriplet[0]]->    colourNeighbour(outgoing[octet[0]]);
//       if(inter&&inter->coloured()) {
// 	if(inter->dataPtr()->iColour()==PDT::Colour3)
// 	  outgoing[antitriplet[0]]->antiColourLine()->addColoured(inter);
// 	else if(inter->dataPtr()->iColour()==PDT::Colour3bar)
// 	  outgoing[    triplet[0]]->    colourLine()->addAntiColoured(inter);
// 	else if(inter->dataPtr()->iColour()==PDT::Colour8) {
// 	  outgoing[antitriplet[0]]->antiColourLine()->addAntiColoured(inter);
// 	  outgoing[    triplet[0]]->    colourLine()->addColoured(inter);
// 	}
//       }
//     }
//     else if(triplet.size()==3) {
//       tColinePtr col[3] = {ColourLine::create(outgoing[0]),
// 			   ColourLine::create(outgoing[1]),
// 			   ColourLine::create(outgoing[2])};
//       col[0]->setSourceNeighbours(col[1],col[2]);
//     }
//     else if(antitriplet.size()==3) {
//       tColinePtr col[3] = {ColourLine::create(outgoing[0],true),
// 			   ColourLine::create(outgoing[1],true),
// 			   ColourLine::create(outgoing[2],true)};
//       col[0]->setSinkNeighbours(col[1],col[2]);
//     }
//     else {
//       string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
// 	+ out[1]->PDGName() + " " + out[2]->PDGName();
//       throw Exception() 
// 	<< "Unknown colour structure in GeneralFourBodyDecayer::"
// 	<< "colourConnections() for singlet decaying particle "
// 	<< mode << Exception::runerror;
//     } 
//   }
//   // colour triplet decaying particle
//   else if( inColour == PDT::Colour3) {
//     if(singlet.size()==2&&triplet.size()==1) {
//       outgoing[triplet[0]]->incomingColour(const_ptr_cast<tPPtr>(&parent));
//       if(inter&&inter->coloured()) 
// 	outgoing[triplet[0]]->colourLine()->addColoured(inter);
//     }
//     else if(antitriplet.size()==1&&triplet.size()==2) {
//       if(colourFlow()==0) {
// 	outgoing[triplet[0]]->incomingColour(const_ptr_cast<tPPtr>(&parent));
// 	outgoing[antitriplet[0]]->colourNeighbour(outgoing[triplet[1]]);
// 	if(inter&&inter->coloured()) {
// 	  switch (inter->dataPtr()->iColour()) {
// 	  case PDT::Colour8:
// 	    inter->incomingColour(const_ptr_cast<tPPtr>(&parent));
// 	    outgoing[triplet[1]]->colourLine()->addAntiColoured(inter);
// 	    break;
// 	  default:
// 	    string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
// 	      + out[1]->PDGName() + " " + out[2]->PDGName();
// 	    throw Exception() << "Unknown colour for intermediate in "
// 			      << "GeneralFourBodyDecayer::"
// 			      << "colourConnections() for "
// 			      << "decaying colour triplet " 
// 			      << mode << Exception::runerror;
// 	  }
// 	}
//       }
//       else {
// 	outgoing[triplet[1]]->incomingColour(const_ptr_cast<tPPtr>(&parent));
// 	outgoing[antitriplet[0]]->colourNeighbour(outgoing[triplet[0]]);
// 	if(inter&&inter->coloured()) {
// 	  switch (inter->dataPtr()->iColour()) {
// 	  case PDT::Colour8:
// 	    inter->incomingColour(const_ptr_cast<tPPtr>(&parent));
// 	    outgoing[triplet[0]]->colourLine()->addAntiColoured(inter);
// 	    break;
// 	  default:
// 	    string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
// 	      + out[1]->PDGName() + " " + out[2]->PDGName();
// 	    throw Exception() << "Unknown colour for intermediate in "
// 			      << "GeneralFourBodyDecayer::"
// 			      << "colourConnections() for "
// 			      << "decaying colour triplet " 
// 			      << mode << Exception::runerror;
// 	  }
// 	}
//       }
//     }
//     else if (singlet.size()==1&&triplet.size()==1&&octet.size()==1) {
//       if(inter) {
// 	if(inter->children()[0]->dataPtr()->iColour()==PDT::Colour8 ||
// 	   inter->children()[1]->dataPtr()->iColour()==PDT::Colour8) {
// 	  inter->incomingColour(const_ptr_cast<tPPtr>(&parent));
// 	  outgoing[octet[0]]->incomingColour(inter);
// 	  outgoing[octet[0]]->colourNeighbour(outgoing[triplet[0]]);
// 	}
// 	else {
// 	  outgoing[octet[0]]->incomingColour(inter);
// 	  outgoing[octet[0]]->colourNeighbour(inter);
// 	  outgoing[triplet[0]]->incomingColour(inter);
// 	}
//       }
//       else {
// 	outgoing[octet[0]]->incomingColour(const_ptr_cast<tPPtr>(&parent));
// 	outgoing[octet[0]]->colourNeighbour(outgoing[triplet[0]]);
//       }
//     }
//     else {
//       string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
// 	+ out[1]->PDGName() + " " + out[2]->PDGName();
//       throw Exception() 
// 	<< "Unknown colour structure in GeneralFourBodyDecayer::"
// 	<< "colourConnections() for triplet decaying particle " 
// 	<< mode << Exception::runerror;
//     }
//   }
//   else if( inColour == PDT::Colour3bar) {
//     if(singlet.size()==2&&antitriplet.size()==1) {
//       outgoing[antitriplet[0]]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
//     }
//     else if(antitriplet.size()==2&&triplet.size()==1) {
//       if(colourFlow()==0) {
// 	outgoing[antitriplet[0]]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
// 	outgoing[triplet[0]]->antiColourNeighbour(outgoing[antitriplet[1]]);
// 	if(inter&&inter->coloured()) {
// 	  switch (inter->dataPtr()->iColour()) {
// 	  case PDT::Colour8:
// 	    inter->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
// 	    outgoing[antitriplet[1]]->antiColourLine()->addAntiColoured(inter);
// 	    break;
// 	  default:
// 	    string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
// 	      + out[1]->PDGName() + " " + out[2]->PDGName();
// 	    throw Exception() << "Unknown colour for intermediate in"
// 			      << " GeneralFourBodyDecayer::"
// 			      << "colourConnections() for "
// 			      << "decaying colour antitriplet " 
// 			      << mode << Exception::runerror;
// 	  }
// 	}
//       }
//       else {
// 	outgoing[antitriplet[1]]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
// 	outgoing[triplet[0]]->antiColourNeighbour(outgoing[antitriplet[0]]);
// 	if(inter&&inter->coloured()) {
// 	  switch (inter->dataPtr()->iColour()) {
// 	  case PDT::Colour8:
// 	    inter->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
// 	    outgoing[antitriplet[0]]->antiColourLine()->addAntiColoured(inter);
// 	    break;
// 	  default:
// 	    string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
// 	      + out[1]->PDGName() + " " + out[2]->PDGName();
// 	    throw Exception() << "Unknown colour for intermediate in "
// 			      << "GeneralFourBodyDecayer::"
// 			      << "colourConnections() for "
// 			      << "decaying colour antitriplet " 
// 			      << mode << Exception::runerror;
// 	  }
// 	}
//       }
//     }
//     else if (singlet.size()==1&&antitriplet.size()==1&&octet.size()==1) {
//       if(inter) {
// 	if(inter->children()[0]->dataPtr()->iColour()==PDT::Colour8 ||
// 	   inter->children()[1]->dataPtr()->iColour()==PDT::Colour8) {
// 	  inter->incomingColour(const_ptr_cast<tPPtr>(&parent));
// 	  outgoing[octet[0]]->incomingAntiColour(inter);
// 	  outgoing[octet[0]]->antiColourNeighbour(outgoing[antitriplet[0]]);
// 	}
// 	else {
// 	  outgoing[octet[0]]->incomingAntiColour(inter);
// 	  outgoing[octet[0]]->antiColourNeighbour(inter);
// 	  outgoing[antitriplet[0]]->incomingAntiColour(inter);
// 	}
//       }
//       else {
// 	outgoing[octet[0]]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
// 	outgoing[octet[0]]->antiColourNeighbour(outgoing[antitriplet[0]]);
//       }
//     }
//     else {
//       string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
// 	+ out[1]->PDGName() + " " + out[2]->PDGName();
//       throw Exception() 
// 	<< "Unknown colour structure in GeneralFourBodyDecayer::"
// 	<< "colourConnections() for anti-triplet decaying particle" 
// 	<< mode << Exception::runerror;
//     }
//   }
//   else if( inColour == PDT::Colour8) {
//     if(triplet.size()==1&&antitriplet.size()==1&&singlet.size()==1) {
//       outgoing[    triplet[0]]->incomingColour    (const_ptr_cast<tPPtr>(&parent));
//       outgoing[antitriplet[0]]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
//       if(inter&&inter->coloured()) {
// 	switch (inter->dataPtr()->iColour()) {
// 	case PDT::Colour3:
// 	  outgoing[triplet[0]]->colourLine()->addColoured(inter);
// 	  break;
// 	case PDT::Colour3bar:
// 	  outgoing[antitriplet[0]]->antiColourLine()->addAntiColoured(inter);
// 	  break;
// 	default:
// 	  string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
// 	    + out[1]->PDGName() + " " + out[2]->PDGName();
// 	  throw Exception() << "Unknown colour for intermediate"
// 			    << " in GeneralFourBodyDecayer::"
// 			    << "colourConnections() for "
// 			    << "decaying colour octet " 
// 			    << mode << Exception::runerror;
// 	}
//       }
//     }
//     else {
//       string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
// 	+ out[1]->PDGName() + " " + out[2]->PDGName();
//       throw Exception() 
// 	<< "Unknown colour structure in GeneralFourBodyDecayer::"
// 	<< "colourConnections() for octet decaying particle" 
// 	<< mode << Exception::runerror;
//     }
//   }
}

bool GeneralFourBodyDecayer::setColourFactors() {
  string name = _incoming->PDGName() + "->";
  vector<int> sng,trip,atrip,oct;
  unsigned int iloc(0);
  for(vector<PDPtr>::const_iterator it = _outgoing.begin();
      it != _outgoing.end();++it) {
    name += (**it).PDGName() + " ";
    if     ((**it).iColour() == PDT::Colour0    ) sng.push_back(iloc) ;
    else if((**it).iColour() == PDT::Colour3    ) trip.push_back(iloc) ;
    else if((**it).iColour() == PDT::Colour3bar ) atrip.push_back(iloc);
    else if((**it).iColour() == PDT::Colour8    ) oct.push_back(iloc) ;
    ++iloc;
  }
  // colour neutral decaying particle
  if     ( _incoming->iColour() == PDT::Colour0) {
    // options are all neutral or triplet/antitriplet+ neutral
    if(sng.size()==4) {
      _nflow = 1;
      _colour        = vector<DVector>(1,DVector(1,1.));
      _colourLargeNC = vector<DVector>(1,DVector(1,1.));
    }
    else if(sng.size()==2&&trip.size()==1&&atrip.size()==1) {
      _nflow = 1;
      _colour         = vector<DVector>(1,DVector(1,3.));
      _colourLargeNC  = vector<DVector>(1,DVector(1,3.));
    }
//     else if(trip.size()==1&&atrip.size()==1&&oct.size()==1) {
//       _nflow = 1;
//       _colour         = vector<DVector>(1,DVector(1,4.));
//       _colourLargeNC  = vector<DVector>(1,DVector(1,4.));
//     }
//     else if( trip.size() == 3 || atrip.size() == 3 ) {
//       _nflow = 1;
//       _colour         = vector<DVector>(1,DVector(1,6.));
//       _colourLargeNC  = vector<DVector>(1,DVector(1,6.));
//       for(unsigned int ix=0;ix<_diagrams.size();++ix) {
// 	tPDPtr inter = _diagrams[ix].intermediate;
// 	if(inter->CC()) inter = inter->CC();
// 	unsigned int io[2]={1,2};
// 	double sign = _diagrams[ix].channelType == TBDiagram::channel13 ? -1. : 1.;
// 	for(unsigned int iy=0;iy<3;++iy) {
// 	  if     (iy==1) io[0]=0;
// 	  else if(iy==2) io[1]=1;
// 	  tPDVector decaylist = _diagrams[ix].vertices.second->search(iy, inter);
// 	  if(decaylist.empty()) continue;
// 	  bool found=false;
// 	  for(unsigned int iz=0;iz<decaylist.size();iz+=3) {	    
// 	    if(decaylist[iz+io[0]]->id()==_diagrams[ix].outgoingPair.first &&
// 	       decaylist[iz+io[1]]->id()==_diagrams[ix].outgoingPair.second) {
// 	      sign *= 1.;
// 	      found = true;
// 	    }
// 	    else if(decaylist[iz+io[0]]->id()==_diagrams[ix].outgoingPair.second &&
// 		    decaylist[iz+io[1]]->id()==_diagrams[ix].outgoingPair.first ) {
// 	      sign *= -1.;
// 	      found = true;
// 	    }
// 	  }
// 	  if(found) {
// 	    if(iy==1) sign *=-1.;
// 	    break;
// 	  }
// 	}
// 	_diagrams[ix].       colourFlow = vector<CFPair>(1,make_pair(1,sign));
// 	_diagrams[ix].largeNcColourFlow = vector<CFPair>(1,make_pair(1,sign));
//       }
//     }
    else {
      generator()->log() << "Unknown colour flow structure for "
			 << "colour neutral decay " << name 
			 << " in getColourFactors(), omitting decay\n";
      return false;
    }
  }
  // colour triplet decaying particle
  else if( _incoming->iColour() == PDT::Colour3) {
    if(sng.size()==3&&trip.size()==1) {
      _nflow = 1;
      _colour        = vector<DVector>(1,DVector(1,1.));
      _colourLargeNC = vector<DVector>(1,DVector(1,1.));
    }
    else if(trip.size()==2&&atrip.size()==1&&sng.size()==1) {
      _nflow = 2;
      _colour.resize(2,DVector(2,0.));
      _colour[0][0] = 3.; _colour[0][1] = 1.;
      _colour[1][0] = 1.; _colour[1][1] = 3.;
      _colourLargeNC.resize(2,DVector(2,0.));
      _colourLargeNC[0][0] = 3.; _colourLargeNC[1][1] = 3.;
      // particle connect to incoming in first flow
      unsigned int ifirst(0);
      for(unsigned int ix=0;ix<_outgoing.size();++ix) {
	if(_outgoing[ix]->iColour()==PDT::Colour3&&ifirst==0) {
	  ifirst = ix+1;
	  break;
	}
      }
      // sort out the contribution of the different diagrams
      // to the colour flows
      for(vector<NBDiagram>::iterator it = _diagrams.begin();
	  it!=_diagrams.end();++it) {
	// first topology
	if(it->vertices.begin()->second.incoming) {
	  tPDPair inter = make_pair(it->vertices.begin() ->second.incoming,
				    (++(it->vertices.begin()))->second.incoming);
	  // one neutral and one triplet
	  if(inter. first->iColour()==PDT::Colour3 &&
	     inter.second->iColour()==PDT::Colour0) {
	    if( it->channelType/1000       == ifirst ||
	       (it->channelType%1000)/100  == ifirst) {
	      it->       colourFlow = vector<CFPair>(1,make_pair(1,1.));
	      it->largeNcColourFlow = vector<CFPair>(1,make_pair(1,1.));
	    }
	    else {
	      it->colourFlow        = vector<CFPair>(1,make_pair(2,1.));
	      it->largeNcColourFlow = vector<CFPair>(1,make_pair(2,1.));
	    }
	  }
	  else if(inter. first->iColour()==PDT::Colour0 &&
		  inter.second->iColour()==PDT::Colour3) {
	    if((it->channelType%100)/10 == ifirst ||
	        it->channelType%10      == ifirst) {
	      it->       colourFlow = vector<CFPair>(1,make_pair(1,1.));
	      it->largeNcColourFlow = vector<CFPair>(1,make_pair(1,1.));
	    }
	    else {
	      it->colourFlow        = vector<CFPair>(1,make_pair(2,1.));
	      it->largeNcColourFlow = vector<CFPair>(1,make_pair(2,1.));
	    }
	  }
	  // one neutral and one octet
	  else if(inter. first->iColour()==PDT::Colour3 &&
		  inter.second->iColour()==PDT::Colour8) {
	    if( it->channelType/1000       == ifirst ||
	       (it->channelType%1000)/100  == ifirst) {
	      vector<CFPair> flow(1,make_pair(2, 0.5  ));
	      it->largeNcColourFlow = flow;
	      flow.push_back(       make_pair(1,-1./6.));
	      it->colourFlow=flow;
	    }
	    else {
	      vector<CFPair> flow(1,make_pair(1, 0.5  ));
	      it->largeNcColourFlow = flow;
	      flow.push_back(       make_pair(2,-1./6.));
	      it->colourFlow=flow;
	    }
	  }
	  else if(inter. first->iColour()==PDT::Colour8 &&
		  inter.second->iColour()==PDT::Colour3) {
	    if((it->channelType%100)/10 == ifirst ||
	        it->channelType%10      == ifirst) {
	      vector<CFPair> flow(1,make_pair(2, 0.5  ));
	      it->largeNcColourFlow = flow;
	      flow.push_back(       make_pair(1,-1./6.));
	      it->colourFlow=flow;
	    }
	    else {
	      vector<CFPair> flow(1,make_pair(1, 0.5  ));
	      it->largeNcColourFlow = flow;
	      flow.push_back(       make_pair(2,-1./6.));
	      it->colourFlow=flow;
	    }
	  }
	  else
	    return false;
	}
	// second topology
	else {
	  tPDPtr inter = (++(it->vertices.begin()))->second.incoming;
	  unsigned int iloc;
	  if(inter->iColour()==PDT::Colour3) {
	    iloc = (it->channelType%1000)/100;
	    const NBVertex & vertex =  (++(it->vertices.begin()))->second;
	    inter = (++vertex.vertices.begin())->second.incoming;
	  }
	  else {
	    iloc = it->channelType/1000;
	  }
	  if(inter->iColour()==PDT::Colour0) {
	    if( iloc == ifirst) {
	      it->       colourFlow = vector<CFPair>(1,make_pair(1,1.));
	      it->largeNcColourFlow = vector<CFPair>(1,make_pair(1,1.));
	    }
	    else {
	      it->colourFlow        = vector<CFPair>(1,make_pair(2,1.));
	      it->largeNcColourFlow = vector<CFPair>(1,make_pair(2,1.));
	    }
	  }
	  else if(inter->iColour()==PDT::Colour8) {
	    if( iloc == ifirst) {
	      vector<CFPair> flow(1,make_pair(2, 0.5  ));
	      it->largeNcColourFlow = flow;
	      flow.push_back(       make_pair(1,-1./6.));
	      it->colourFlow=flow;
	    }
	    else {
	      vector<CFPair> flow(1,make_pair(1, 0.5  ));
	      it->largeNcColourFlow = flow;
	      flow.push_back(       make_pair(2,-1./6.));
	      it->colourFlow=flow;
	    }
	  }
	  else
	    return false;
	}
      }
    }
    else {
      generator()->log() << "Unknown colour structure for "
			 << "triplet decay in "
			 << "FourBodyDecayConstructor::getColourFactors()"
			 << " for " << name << " omitting decay\n";
      return false;
    }
  }
  // colour antitriplet decaying particle
  else if( _incoming->iColour() == PDT::Colour3bar) {
    if(sng.size()==3&&atrip.size()==1) {
      _nflow = 1;
      _colour        = vector<DVector>(1,DVector(1,1.));
      _colourLargeNC = vector<DVector>(1,DVector(1,1.));
    }
    else if(atrip.size()==2&&trip.size()==1&&sng.size()==1) {
      _nflow = 2;
      _colour.resize(2,DVector(2,0.));
      _colour[0][0] = 3.; _colour[0][1] = 1.;
      _colour[1][0] = 1.; _colour[1][1] = 3.;
      _colourLargeNC.resize(2,DVector(2,0.));
      // particle connect to incoming in first flow
      unsigned int ifirst(0);
      for(unsigned int ix=0;ix<_outgoing.size();++ix) {
	if(_outgoing[ix]->iColour()==PDT::Colour3&&ifirst==0) {
	  ifirst = ix+1;
	  break;
	}
      }
      // sort out the contribution of the different diagrams
      // to the colour flows
      for(vector<NBDiagram>::iterator it = _diagrams.begin();
	  it!=_diagrams.end();++it) {
	// first topology
	if(it->vertices.begin()->second.incoming) {
	  tPDPair inter = make_pair(it->vertices.begin() ->second.incoming,
				    (++(it->vertices.begin()))->second.incoming);
	  // one neutral and one triplet
	  if(inter. first->iColour()==PDT::Colour3bar &&
	     inter.second->iColour()==PDT::Colour0) {
	    if( it->channelType/1000       == ifirst ||
	       (it->channelType%1000)/100  == ifirst) {
	      it->       colourFlow = vector<CFPair>(1,make_pair(1,1.));
	      it->largeNcColourFlow = vector<CFPair>(1,make_pair(1,1.));
	    }
	    else {
	      it->colourFlow        = vector<CFPair>(1,make_pair(2,1.));
	      it->largeNcColourFlow = vector<CFPair>(1,make_pair(2,1.));
	    }
	  }
	  else if(inter. first->iColour()==PDT::Colour0 &&
		  inter.second->iColour()==PDT::Colour3bar) {
	    if((it->channelType%100)/10 == ifirst ||
	        it->channelType%10      == ifirst) {
	      it->       colourFlow = vector<CFPair>(1,make_pair(1,1.));
	      it->largeNcColourFlow = vector<CFPair>(1,make_pair(1,1.));
	    }
	    else {
	      it->colourFlow        = vector<CFPair>(1,make_pair(2,1.));
	      it->largeNcColourFlow = vector<CFPair>(1,make_pair(2,1.));
	    }
	  }
	  // one neutral and one octet
	  else if(inter. first->iColour()==PDT::Colour3bar &&
		  inter.second->iColour()==PDT::Colour8) {
	    if( it->channelType/1000       == ifirst ||
	       (it->channelType%1000)/100  == ifirst) {
	      vector<CFPair> flow(1,make_pair(2, 0.5  ));
	      it->largeNcColourFlow = flow;
	      flow.push_back(       make_pair(1,-1./6.));
	      it->colourFlow=flow;
	    }
	    else {
	      vector<CFPair> flow(1,make_pair(1, 0.5  ));
	      it->largeNcColourFlow = flow;
	      flow.push_back(       make_pair(2,-1./6.));
	      it->colourFlow=flow;
	    }
	  }
	  else if(inter. first->iColour()==PDT::Colour8 &&
		  inter.second->iColour()==PDT::Colour3bar) {
	    if((it->channelType%100)/10 == ifirst ||
	        it->channelType%10      == ifirst) {
	      vector<CFPair> flow(1,make_pair(2, 0.5  ));
	      it->largeNcColourFlow = flow;
	      flow.push_back(       make_pair(1,-1./6.));
	      it->colourFlow=flow;
	    }
	    else {
	      vector<CFPair> flow(1,make_pair(1, 0.5  ));
	      it->largeNcColourFlow = flow;
	      flow.push_back(       make_pair(2,-1./6.));
	      it->colourFlow=flow;
	    }
	  }
	  else
	    return false;
	}
	// second topology
	else {
	  tPDPtr inter = (++(it->vertices.begin()))->second.incoming;
	  unsigned int iloc;
	  if(inter->iColour()==PDT::Colour3bar) {
	    iloc = (it->channelType%1000)/100;
	    const NBVertex & vertex =  (++(it->vertices.begin()))->second;
	    inter = (++vertex.vertices.begin())->second.incoming;
	  }
	  else {
	    iloc = it->channelType/1000;
	  }
	  if(inter->iColour()==PDT::Colour0) {
	    if( iloc == ifirst) {
	      it->       colourFlow = vector<CFPair>(1,make_pair(1,1.));
	      it->largeNcColourFlow = vector<CFPair>(1,make_pair(1,1.));
	    }
	    else {
	      it->colourFlow        = vector<CFPair>(1,make_pair(2,1.));
	      it->largeNcColourFlow = vector<CFPair>(1,make_pair(2,1.));
	    }
	  }
	  else if(inter->iColour()==PDT::Colour8) {
	    if( iloc == ifirst) {
	      vector<CFPair> flow(1,make_pair(2, 0.5  ));
	      it->largeNcColourFlow = flow;
	      flow.push_back(       make_pair(1,-1./6.));
	      it->colourFlow=flow;
	    }
	    else {
	      vector<CFPair> flow(1,make_pair(1, 0.5  ));
	      it->largeNcColourFlow = flow;
	      flow.push_back(       make_pair(2,-1./6.));
	      it->colourFlow=flow;
	    }
	  }
	  else
	    return false;
	}
      }
    }
    else {
      generator()->log() << "Unknown colour antitriplet decay in "
			 << "FourBodyDecayConstructor::getColourFactors()"
			 << " for " << name << " omitting decay\n";
      return false;
    }
  }
  else if( _incoming->iColour() == PDT::Colour8) {
    // triplet antitriplet
    if(trip.size() == 1 && atrip.size() == 1 && sng.size() == 2) {
      _nflow = 1;
      _colour        = vector<DVector>(1,DVector(1,0.5));
      _colourLargeNC = vector<DVector>(1,DVector(1,0.5));
    }
    else {
      generator()->log() << "Unknown colour octet decay in "
			 << "FourBodyDecayConstructor::getColourFactors()"
			 << " for " << name << " omitting decay\n";
      return false;
    }
  }
  if( Debug::level > 1 ) {
    generator()->log() << "Mode: " << name << " has colour factors\n";
    for(unsigned int ix=0;ix<_nflow;++ix) {
      for(unsigned int iy=0;iy<_nflow;++iy) {
	generator()->log() << _colour[ix][iy] << " ";
      }
      generator()->log() << "\n";
    }
  }
  return true;
}




















