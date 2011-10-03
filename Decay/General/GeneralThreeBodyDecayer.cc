// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralThreeBodyDecayer class.
//

#include "GeneralThreeBodyDecayer.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/PDT/ThreeBodyAllOnCalculator.h"

using namespace Herwig;

/**
 *  A struct to order the particles in the same way as in the DecayMode's
 */
struct ParticleOrdering {
  bool operator()(PDPtr p1, PDPtr p2) {
    return abs(p1->id()) > abs(p2->id()) ||
      ( abs(p1->id()) == abs(p2->id()) && p1->id() > p2->id() ) ||
      ( p1->id() == p2->id() && p1->fullName() > p2->fullName() );
  }
};

/**
 * A set of ParticleData objects ordered as for the DecayMode's
 */
typedef multiset<PDPtr,ParticleOrdering> OrderedParticles;

void GeneralThreeBodyDecayer::persistentOutput(PersistentOStream & os) const {
  os << _incoming << _outgoing << _diagrams << _diagmap << _colour << _colourLargeNC
     << _nflow << _widthopt << _reftag << _reftagcc << _intOpt << _relerr;
}

void GeneralThreeBodyDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incoming >> _outgoing >> _diagrams >> _diagmap >> _colour >> _colourLargeNC
     >> _nflow >> _widthopt >> _reftag >> _reftagcc >> _intOpt >> _relerr;
}

AbstractClassDescription<GeneralThreeBodyDecayer> 
GeneralThreeBodyDecayer::initGeneralThreeBodyDecayer;
// Definition of the static class description member.

void GeneralThreeBodyDecayer::Init() {

  static ClassDocumentation<GeneralThreeBodyDecayer> documentation
    ("The GeneralThreeBodyDecayer class is the base class for the implementation of"
     " all three body decays based on spin structures in Herwig++.");

  static Switch<GeneralThreeBodyDecayer,unsigned int> interfaceWidthOption
    ("WidthOption",
     "Option for the treatment of the widths of the intermediates",
     &GeneralThreeBodyDecayer::_widthopt, 1, false, false);
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

  static Switch<GeneralThreeBodyDecayer,unsigned int> interfacePartialWidthIntegration
    ("PartialWidthIntegration",
     "Switch to control the partial width integration",
     &GeneralThreeBodyDecayer::_intOpt, 0, false, false);
  static SwitchOption interfacePartialWidthIntegrationAllPoles
    (interfacePartialWidthIntegration,
     "AllPoles",
     "Include all potential poles",
     0);
  static SwitchOption interfacePartialWidthIntegrationShallowestPole
    (interfacePartialWidthIntegration,
     "ShallowestPole",
     "Only include the shallowest pole",
     1);

  static Parameter<GeneralThreeBodyDecayer,double> interfaceRelativeError
    ("RelativeError",
     "The relative error for the GQ integration of the partial width",
     &GeneralThreeBodyDecayer::_relerr, 1e-2, 1e-10, 1.,
     false, false, Interface::limited);

}

ParticleVector GeneralThreeBodyDecayer::decay(const Particle & parent,
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

int  GeneralThreeBodyDecayer::
modeNumber(bool & cc, tcPDPtr in, const tPDVector & outin) const {
  assert( !_reftag.empty() && !_reftagcc.empty() );
  // check number of outgoing particles
  if( outin.size() != 3 || abs(in->id()) != _incoming->id() ) return -1;
  OrderedParticles testmode(outin.begin(), outin.end());
  OrderedParticles::const_iterator dit = testmode.begin();
  string testtag(in->name() + "->");
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

void GeneralThreeBodyDecayer::setDecayInfo(PDPtr incoming,
					   vector<PDPtr> outgoing,
					   const vector<TBDiagram> & process,
					   const vector<DVector> & factors,
					   const vector<DVector> & Ncfactors,
					   const unsigned int ncf) {
  // set the member variables from the info supplied
  _incoming        = incoming;
  _outgoing        = outgoing;
  _diagrams        = process;
  _colour          = factors;
  _colourLargeNC   = Ncfactors;
  _nflow           = ncf;
  assert( _outgoing.size() == 3 );
  // Construct reference tags for testing in modeNumber function
  OrderedParticles refmode(_outgoing.begin(), _outgoing.end());
  OrderedParticles::const_iterator dit = refmode.begin();
  _reftag = _incoming->name() + "->";
  for( unsigned int i = 1; dit != refmode.end(); ++dit, ++i) {
    _reftag += (**dit).name();
    if( i != 3 )  _reftag += string(",");
  }
  //CC-mode
  refmode.clear();
  _reftagcc = _incoming->CC() ? _incoming->CC()->name() : 
    _incoming->name();
  _reftagcc += "->";
  for( unsigned int i = 0;  i < 3; ++i ) {
    if( _outgoing[i]->CC() ) refmode.insert( _outgoing[i]->CC() );
    else refmode.insert( _outgoing[i] );
  }
  dit = refmode.begin();
  for( unsigned int i = 1; dit != refmode.end(); ++dit , ++i) {
    _reftagcc += (**dit).name();
    if( i != 3 ) _reftagcc += string(",");
  }
}

void GeneralThreeBodyDecayer::doinit() {
  DecayIntegrator::doinit();
  // create the phase space integrator
  tPDVector extpart(1,_incoming);
  extpart.insert(extpart.end(),_outgoing.begin(),_outgoing.end());
  // create the integration channels for the decay
  DecayPhaseSpaceModePtr mode(new_ptr(DecayPhaseSpaceMode(extpart,this,true)));
  DecayPhaseSpaceChannelPtr newchannel;
  // create the phase-space channels for the integration
  unsigned int nmode(0);
  for(unsigned int ix=0;ix<_diagrams.size();++ix) {
    if(_diagrams[ix].channelType==TBDiagram::fourPoint||
       _diagrams[ix].channelType==TBDiagram::UNDEFINED) continue;
    // create the new channel
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    if(_diagrams[ix].channelType==TBDiagram::channel23) {
      newchannel->addIntermediate(extpart[0],0,0.0,-1,1);
      newchannel->addIntermediate(_diagrams[ix].intermediate,0,0.0, 2,3);
    }
    else if(_diagrams[ix].channelType==TBDiagram::channel13) {
      newchannel->addIntermediate(extpart[0],0,0.0,-1,2);
      newchannel->addIntermediate(_diagrams[ix].intermediate,0,0.0, 1,3);
    }
    else if(_diagrams[ix].channelType==TBDiagram::channel12) {
      newchannel->addIntermediate(extpart[0],0,0.0,-1,3);
      newchannel->addIntermediate(_diagrams[ix].intermediate,0,0.0, 1,2);
    }
    _diagmap.push_back(ix);
    mode->addChannel(newchannel);
    ++nmode;
  }
  if(nmode==0) {
    string mode = extpart[0]->PDGName() + "->";
    for(unsigned int ix=1;ix<extpart.size();++ix) mode += extpart[ix]->PDGName() + " ";
    throw Exception() << "No decay channels in GeneralThreeBodyDecayer::"
		      << "doinit() for " << mode << "\n" << Exception::runerror;
  }
  // add the mode
  vector<double> wgt(nmode,1./double(nmode));
  addMode(mode,1.,wgt);
}

double GeneralThreeBodyDecayer::
threeBodyMatrixElement(const int imode,  const Energy2 q2,
		       const Energy2 s3, const Energy2 s2, 
		       const Energy2 s1, const Energy  m1, 
		       const Energy  m2, const Energy  m3) const {
  // calculate the momenta of the outgoing particles
  Energy m0=sqrt(q2);
  // energies
  Energy eout[3] = {0.5*(q2+sqr(m1)-s1)/m0,
		    0.5*(q2+sqr(m2)-s2)/m0,
		    0.5*(q2+sqr(m3)-s3)/m0};
  // magnitudes of the momenta
  Energy pout[3] = {sqrt(sqr(eout[0])-sqr(m1)),
		    sqrt(sqr(eout[1])-sqr(m2)),
		    sqrt(sqr(eout[2])-sqr(m3))};
  double cos2 = 0.5*(sqr(pout[0])+sqr(pout[1])-sqr(pout[2]))/pout[0]/pout[1];
  double cos3 = 0.5*(sqr(pout[0])-sqr(pout[1])+sqr(pout[2]))/pout[0]/pout[2];
  double sin2 = sqrt(1.-sqr(cos2)), sin3 = sqrt(1.-sqr(cos3));
  Lorentz5Momentum out[3]=
    {Lorentz5Momentum(      ZERO   , ZERO ,  pout[0]      , eout[0] , m1),
     Lorentz5Momentum(  pout[1]*sin2 , ZERO , -pout[1]*cos2 , eout[1] , m2),
     Lorentz5Momentum( -pout[2]*sin3 , ZERO , -pout[2]*cos3 , eout[2] , m3)};
  // create the incoming
  PPtr inpart=mode(imode)->externalParticles(0)->
    produceParticle(Lorentz5Momentum(sqrt(q2)));
  // and outgoing particles
  ParticleVector decay;
  for(unsigned int ix=1;ix<4;++ix)
    decay.push_back(mode(imode)->externalParticles(ix)->produceParticle(out[ix-1]));
  // return the matrix element
  return me2(-1,*inpart,decay,Initialize);
}

double GeneralThreeBodyDecayer::brat(const DecayMode &, const Particle & p,
				     double oldbrat) const {
  ParticleVector children = p.children();
  if( children.size() != 3 || !p.data().widthGenerator() ) 
    return oldbrat;
  // partial width for this mode
  Energy scale = p.mass();
  Energy pwidth = 
    partialWidth( make_pair(p.dataPtr(), scale),
		  make_pair(children[0]->dataPtr(), children[0]->mass()),
		  make_pair(children[1]->dataPtr(), children[1]->mass()),
		  make_pair(children[2]->dataPtr(), children[2]->mass()) );
  Energy width = p.data().widthGenerator()->width(p.data(), scale);
  return pwidth/width;
}

Energy GeneralThreeBodyDecayer::partialWidth(PMPair inpart, PMPair outa, 
					     PMPair outb, PMPair outc) const {
  if(inpart.second<outa.second+outb.second+outc.second) return ZERO;
  // create the object to calculate the width if it doesn't all ready exist
  if(!_widthcalc) {
    string tag = _incoming->name() + "->";
    tag += _outgoing[0]->name() + "," + _outgoing[1]->name() + ","
      + _outgoing[2]->name() + ";";
    DMPtr dm = generator()->findDecayMode(tag);
    _widthcalc = threeBodyMEIntegrator(*dm);
  }
  return _widthcalc->partialWidth(sqr(inpart.second));
}

void GeneralThreeBodyDecayer::
colourConnections(const Particle & parent,
		  const ParticleVector & out) const {
  // first extract the outgoing particles and intermediate
  PPtr inter;
  ParticleVector outgoing;
  if(!generateIntermediates()) {
    outgoing=out;
  }
  else {
    // find the diagram
    unsigned int idiag = diagramMap()[mode(imode())->selectedChannel()];
    PPtr child;
    for(unsigned int ix=0;ix<out.size();++ix) {
      if(out[ix]->children().empty()) child = out[ix];
      else                            inter = out[ix];
    }
    outgoing.resize(3);
    switch(_diagrams[idiag].channelType) {
    case TBDiagram::channel23:
      outgoing[0] = child;
      outgoing[1] = inter->children()[0];
      outgoing[2] = inter->children()[1];
      break;
    case TBDiagram::channel13:
      outgoing[0] = inter->children()[0];
      outgoing[1] = child;
      outgoing[2] = inter->children()[1];
      break;
    case TBDiagram::channel12:
      outgoing[0] = inter->children()[0];
      outgoing[1] = inter->children()[1];
      outgoing[2] = child;
      break;
    default:
      throw Exception() << "unknown diagram type in GeneralThreeBodyDecayer::"
			<< "colourConnections()" << Exception::runerror;
    }
  }
  // extract colour of the incoming and outgoing particles
  PDT::Colour inColour(parent.data().iColour());
  vector<PDT::Colour> outColour;
  vector<int> singlet,octet,triplet,antitriplet;
  for(unsigned int ix=0;ix<outgoing.size();++ix) {
    outColour.push_back(outgoing[ix]->data().iColour());
    switch(outColour.back()) {
    case PDT::Colour0   :     
      singlet.push_back(ix);
      break;
    case PDT::Colour3   :     
      triplet.push_back(ix);
      break;
    case PDT::Colour3bar: 
      antitriplet.push_back(ix);
      break;
    case PDT::Colour8   :     
      octet.push_back(ix);
      break;
    default:
      throw Exception() << "Unknown colour for particle in GeneralThreeBodyDecayer::"
			<< "colourConnections()" << Exception::runerror;
    }
  }
  // colour neutral decaying particle
  if     ( inColour == PDT::Colour0) {
    // options are all neutral or triplet/antitriplet+ neutral
    if(singlet.size()==3) return;
    else if(singlet.size()==1&&triplet.size()==1&&antitriplet.size()==1) {
      outgoing[triplet[0]]->antiColourNeighbour(outgoing[antitriplet[0]]);
      // add intermediate if needed
      if(inter&&inter->coloured()) {
	if(inter->dataPtr()->iColour()==PDT::Colour3)
	  outgoing[triplet[0]]->colourLine()->addColoured(inter);
	else if(inter->dataPtr()->iColour()==PDT::Colour3bar)
	  outgoing[triplet[0]]->colourLine()->addAntiColoured(inter);
      }
    }
    else if(octet.size()==1&&triplet.size()==1&&antitriplet.size()==1) {
      outgoing[    triplet[0]]->antiColourNeighbour(outgoing[octet[0]]);
      outgoing[antitriplet[0]]->    colourNeighbour(outgoing[octet[0]]);
      if(inter&&inter->coloured()) {
	if(inter->dataPtr()->iColour()==PDT::Colour3)
	  outgoing[antitriplet[0]]->antiColourLine()->addColoured(inter);
	else if(inter->dataPtr()->iColour()==PDT::Colour3bar)
	  outgoing[    triplet[0]]->    colourLine()->addAntiColoured(inter);
	else if(inter->dataPtr()->iColour()==PDT::Colour8) {
	  outgoing[antitriplet[0]]->antiColourLine()->addAntiColoured(inter);
	  outgoing[    triplet[0]]->    colourLine()->addColoured(inter);
	}
      }
    }
    else {
      string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
	+ out[1]->PDGName() + " " + out[2]->PDGName();
      throw Exception() 
	<< "Unknown colour structure in GeneralThreeBodyDecayer::"
	<< "colourConnections() for singlet decaying particle "
	<< mode << Exception::runerror;
    } 
  }
  // colour triplet decaying particle
  else if( inColour == PDT::Colour3) {
    if(singlet.size()==2&&triplet.size()==1) {
      outgoing[triplet[0]]->incomingColour(const_ptr_cast<tPPtr>(&parent));
      if(inter&&inter->coloured()) 
	outgoing[triplet[0]]->colourLine()->addColoured(inter);
    }
    else if(antitriplet.size()==1&&triplet.size()==2) {
      if(colourFlow()==0) {
	outgoing[triplet[0]]->incomingColour(const_ptr_cast<tPPtr>(&parent));
	outgoing[antitriplet[0]]->colourNeighbour(outgoing[triplet[1]]);
	if(inter&&inter->coloured()) {
	  switch (inter->dataPtr()->iColour()) {
	  case PDT::Colour8:
	    inter->incomingColour(const_ptr_cast<tPPtr>(&parent));
	    outgoing[triplet[1]]->colourLine()->addAntiColoured(inter);
	    break;
	  default:
	    string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
	      + out[1]->PDGName() + " " + out[2]->PDGName();
	    throw Exception() << "Unknown colour for intermediate in "
			      << "GeneralThreeBodyDecayer::"
			      << "colourConnections() for "
			      << "decaying colour triplet " 
			      << mode << Exception::runerror;
	  }
	}
      }
      else {
	outgoing[triplet[1]]->incomingColour(const_ptr_cast<tPPtr>(&parent));
	outgoing[antitriplet[0]]->colourNeighbour(outgoing[triplet[0]]);
	if(inter&&inter->coloured()) {
	  switch (inter->dataPtr()->iColour()) {
	  case PDT::Colour8:
	    inter->incomingColour(const_ptr_cast<tPPtr>(&parent));
	    outgoing[triplet[0]]->colourLine()->addAntiColoured(inter);
	    break;
	  default:
	    string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
	      + out[1]->PDGName() + " " + out[2]->PDGName();
	    throw Exception() << "Unknown colour for intermediate in "
			      << "GeneralThreeBodyDecayer::"
			      << "colourConnections() for "
			      << "decaying colour triplet " 
			      << mode << Exception::runerror;
	  }
	}
      }
    }
    else if (singlet.size()==1&&triplet.size()==1&&octet.size()==1) {
      if(inter) {
	if(inter->children()[0]->dataPtr()->iColour()==PDT::Colour8 ||
	   inter->children()[1]->dataPtr()->iColour()==PDT::Colour8) {
	  inter->incomingColour(const_ptr_cast<tPPtr>(&parent));
	  outgoing[octet[0]]->incomingColour(inter);
	  outgoing[octet[0]]->colourNeighbour(outgoing[triplet[0]]);
	}
	else {
	  outgoing[octet[0]]->incomingColour(inter);
	  outgoing[octet[0]]->colourNeighbour(inter);
	  outgoing[triplet[0]]->incomingColour(inter);
	}
      }
      else {
	outgoing[octet[0]]->incomingColour(const_ptr_cast<tPPtr>(&parent));
	outgoing[octet[0]]->colourNeighbour(outgoing[triplet[0]]);
      }
    }
    else {
      string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
	+ out[1]->PDGName() + " " + out[2]->PDGName();
      throw Exception() 
	<< "Unknown colour structure in GeneralThreeBodyDecayer::"
	<< "colourConnections() for triplet decaying particle " 
	<< mode << Exception::runerror;
    }
  }
  else if( inColour == PDT::Colour3bar) {
    if(singlet.size()==2&&antitriplet.size()==1) {
      outgoing[antitriplet[0]]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
    }
    else if(antitriplet.size()==2&&triplet.size()==1) {
      if(colourFlow()==0) {
	outgoing[antitriplet[0]]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
	outgoing[triplet[0]]->antiColourNeighbour(outgoing[antitriplet[1]]);
	if(inter&&inter->coloured()) {
	  switch (inter->dataPtr()->iColour()) {
	  case PDT::Colour8:
	    inter->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
	    outgoing[antitriplet[1]]->antiColourLine()->addAntiColoured(inter);
	    break;
	  default:
	    string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
	      + out[1]->PDGName() + " " + out[2]->PDGName();
	    throw Exception() << "Unknown colour for intermediate in"
			      << " GeneralThreeBodyDecayer::"
			      << "colourConnections() for "
			      << "decaying colour antitriplet " 
			      << mode << Exception::runerror;
	  }
	}
      }
      else {
	outgoing[antitriplet[1]]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
	outgoing[triplet[0]]->antiColourNeighbour(outgoing[antitriplet[0]]);
	if(inter&&inter->coloured()) {
	  switch (inter->dataPtr()->iColour()) {
	  case PDT::Colour8:
	    inter->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
	    outgoing[antitriplet[0]]->antiColourLine()->addAntiColoured(inter);
	    break;
	  default:
	    string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
	      + out[1]->PDGName() + " " + out[2]->PDGName();
	    throw Exception() << "Unknown colour for intermediate in "
			      << "GeneralThreeBodyDecayer::"
			      << "colourConnections() for "
			      << "decaying colour antitriplet " 
			      << mode << Exception::runerror;
	  }
	}
      }
    }
    else if (singlet.size()==1&&antitriplet.size()==1&&octet.size()==1) {
      if(inter) {
	if(inter->children()[0]->dataPtr()->iColour()==PDT::Colour8 ||
	   inter->children()[1]->dataPtr()->iColour()==PDT::Colour8) {
	  inter->incomingColour(const_ptr_cast<tPPtr>(&parent));
	  outgoing[octet[0]]->incomingAntiColour(inter);
	  outgoing[octet[0]]->antiColourNeighbour(outgoing[antitriplet[0]]);
	}
	else {
	  outgoing[octet[0]]->incomingAntiColour(inter);
	  outgoing[octet[0]]->antiColourNeighbour(inter);
	  outgoing[antitriplet[0]]->incomingAntiColour(inter);
	}
      }
      else {
	outgoing[octet[0]]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
	outgoing[octet[0]]->antiColourNeighbour(outgoing[antitriplet[0]]);
      }
    }
    else {
      string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
	+ out[1]->PDGName() + " " + out[2]->PDGName();
      throw Exception() 
	<< "Unknown colour structure in GeneralThreeBodyDecayer::"
	<< "colourConnections() for anti-triplet decaying particle" 
	<< mode << Exception::runerror;
    }
  }
  else if( inColour == PDT::Colour8) {
    if(triplet.size()==1&&antitriplet.size()==1&&singlet.size()==1) {
      outgoing[    triplet[0]]->incomingColour    (const_ptr_cast<tPPtr>(&parent));
      outgoing[antitriplet[0]]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
      if(inter&&inter->coloured()) {
	switch (inter->dataPtr()->iColour()) {
	case PDT::Colour3:
	  outgoing[triplet[0]]->colourLine()->addColoured(inter);
	  break;
	case PDT::Colour3bar:
	  outgoing[antitriplet[0]]->antiColourLine()->addAntiColoured(inter);
	  break;
	default:
	  string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
	    + out[1]->PDGName() + " " + out[2]->PDGName();
	  throw Exception() << "Unknown colour for intermediate"
			    << " in GeneralThreeBodyDecayer::"
			    << "colourConnections() for "
			    << "decaying colour octet " 
			    << mode << Exception::runerror;
	}
      }
    }
    else {
      string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
	+ out[1]->PDGName() + " " + out[2]->PDGName();
      throw Exception() 
	<< "Unknown colour structure in GeneralThreeBodyDecayer::"
	<< "colourConnections() for octet decaying particle" 
	<< mode << Exception::runerror;
    }
  }
}

void GeneralThreeBodyDecayer::
constructIntegratorChannels(vector<int> & intype, vector<Energy> & inmass,
			    vector<Energy> & inwidth, vector<double> & inpow,
			    vector<double> & inweights) const {
  Energy min = incoming()->mass();
  int nchannel(0);
  pair<int,Energy> imin[4]={make_pair(-1,-1.*GeV),make_pair(-1,-1.*GeV),
			    make_pair(-1,-1.*GeV),make_pair(-1,-1.*GeV)};
  Energy absmin = -1e20*GeV;
  int minType   = -1;
  for(unsigned int iy=0;iy<_diagmap.size();++iy) {
    unsigned int ix=_diagmap[iy];
    if(getProcessInfo()[ix].channelType==TBDiagram::fourPoint) continue;
    Energy dm1(min-getProcessInfo()[ix].intermediate->mass());
    Energy dm2(getProcessInfo()[ix].intermediate->mass());
    int itype(0);
    if     (getProcessInfo()[ix].channelType==TBDiagram::channel23) {
      dm1 -= outgoing()[0]->mass();
      dm2 -= outgoing()[1]->mass()+outgoing()[2]->mass();
      itype = 3;
    }
    else if(getProcessInfo()[ix].channelType==TBDiagram::channel13) {
      dm1 -= outgoing()[1]->mass();
      dm2 -= outgoing()[0]->mass()+outgoing()[2]->mass();
      itype = 2;
    }
    else if(getProcessInfo()[ix].channelType==TBDiagram::channel12) {
      dm1 -= outgoing()[2]->mass();
      dm2 -= outgoing()[0]->mass()+outgoing()[1]->mass();
      itype = 1;
    }
    if(dm1<ZERO||dm2<ZERO) {
      if (imin[itype].first < 0  ||
	  (dm1<ZERO && imin[itype].second < dm1)  ) {
	imin[itype] = make_pair(ix,dm1);
	if(dm1<ZERO&&absmin<dm1) {
	  absmin = dm1;
	  minType = itype;
	}
      }
      continue;
    }
    if(getProcessInfo()[ix].intermediate->id()!=ParticleID::gamma) {
      intype.push_back(itype);
      inpow.push_back(0.);
      inmass.push_back(getProcessInfo()[ix].intermediate->mass());
      inwidth.push_back(getProcessInfo()[ix].intermediate->width());
      ++nchannel;
    }
    else if(getProcessInfo()[ix].intermediate->id()==ParticleID::gamma) {
      intype.push_back(itype);
      inpow.push_back(-2.);
      inmass.push_back(-1.*GeV);
      inwidth.push_back(-1.*GeV);
      ++nchannel;
    }
  }
  // physical poles, use them and return
  if(nchannel>0) {
    inweights = vector<double>(nchannel,1./double(nchannel));
    return;
  }
  // use shallowest pole
  else if(_intOpt==1&&minType>0) {
    intype.push_back(minType);
    if(getProcessInfo()[imin[minType].first].intermediate->id()!=ParticleID::gamma) {
      inpow.push_back(0.);
      inmass.push_back(getProcessInfo()[imin[minType].first].intermediate->mass());
      inwidth.push_back(getProcessInfo()[imin[minType].first].intermediate->width());
    }
    else {
      inpow.push_back(-2.);
      inmass.push_back(-1.*GeV);
      inwidth.push_back(-1.*GeV);
    }
    inweights = vector<double>(1,1.);
    return;
  }
  for(unsigned int ix=1;ix<4;++ix) {
    if(imin[ix].first>=0) {
      intype.push_back(ix);
      if(getProcessInfo()[imin[ix].first].intermediate->id()!=ParticleID::gamma) {
	inpow.push_back(0.);
	inmass.push_back(getProcessInfo()[imin[ix].first].intermediate->mass());
	inwidth.push_back(getProcessInfo()[imin[ix].first].intermediate->width());
      }
      else {
	inpow.push_back(-2.);
	inmass.push_back(-1.*GeV);
	inwidth.push_back(-1.*GeV);
      }
      ++nchannel;
    }
  }
  inweights = vector<double>(nchannel,1./double(nchannel));
}
