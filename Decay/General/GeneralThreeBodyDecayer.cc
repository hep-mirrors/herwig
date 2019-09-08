// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralThreeBodyDecayer class.
//

#include "GeneralThreeBodyDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"

using namespace Herwig;

void GeneralThreeBodyDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << diagrams_ << diagmap_ << colour_ << colourLargeNC_
     << nflow_ << widthOpt_ << refTag_ << refTagCC_ << intOpt_ << relerr_;
}

void GeneralThreeBodyDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> diagrams_ >> diagmap_ >> colour_ >> colourLargeNC_
     >> nflow_ >> widthOpt_ >> refTag_ >> refTagCC_ >> intOpt_ >> relerr_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<GeneralThreeBodyDecayer,DecayIntegrator>
describeHerwigGeneralThreeBodyDecayer("Herwig::GeneralThreeBodyDecayer", "Herwig.so");

void GeneralThreeBodyDecayer::Init() {

  static ClassDocumentation<GeneralThreeBodyDecayer> documentation
    ("The GeneralThreeBodyDecayer class is the base class for the implementation of"
     " all three body decays based on spin structures in Herwig.");

  static Switch<GeneralThreeBodyDecayer,unsigned int> interfaceWidthOption
    ("WidthOption",
     "Option for the treatment of the widths of the intermediates",
     &GeneralThreeBodyDecayer::widthOpt_, 1, false, false);
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
     &GeneralThreeBodyDecayer::intOpt_, 0, false, false);
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
     &GeneralThreeBodyDecayer::relerr_, 1e-2, 1e-10, 1.,
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
  assert( !refTag_.empty() && !refTagCC_.empty() );
  // check number of outgoing particles
  if( outin.size() != 3 || abs(in->id()) != abs(incoming_->id()) ) return -1;
  OrderedParticles testmode(outin.begin(), outin.end());
  OrderedParticles::const_iterator dit = testmode.begin();
  string testtag(in->name() + "->");
  for( unsigned int i = 1; dit != testmode.end(); ++dit, ++i) {
    testtag += (**dit).name();
    if( i != 3 ) testtag += string(",");
  }
  if( testtag == refTag_ ) {
    cc = false;
    return 0;
  }
  else if ( testtag == refTagCC_ ) {
    cc = true;
    return 0;
  }
  else return -1;
}

bool GeneralThreeBodyDecayer::setDecayInfo(PDPtr incoming,
					   vector<PDPtr> outgoing,
					   const vector<TBDiagram> & process,
					   double symfac) {
  // set the member variables from the info supplied
  incoming_        = incoming;
  outgoing_        = outgoing;
  diagrams_        = process;
  assert( outgoing_.size() == 3 );
  // Construct reference tags for testing in modeNumber function
  OrderedParticles refmode(outgoing_.begin(), outgoing_.end());
  OrderedParticles::const_iterator dit = refmode.begin();
  refTag_ = incoming_->name() + "->";
  for( unsigned int i = 1; dit != refmode.end(); ++dit, ++i) {
    refTag_ += (**dit).name();
    if( i != 3 )  refTag_ += string(",");
  }
  //CC-mode
  refmode.clear();
  refTagCC_ = incoming_->CC() ? incoming_->CC()->name() : 
    incoming_->name();
  refTagCC_ += "->";
  for( unsigned int i = 0;  i < 3; ++i ) {
    if( outgoing_[i]->CC() ) refmode.insert( outgoing_[i]->CC() );
    else refmode.insert( outgoing_[i] );
  }
  dit = refmode.begin();
  for( unsigned int i = 1; dit != refmode.end(); ++dit , ++i) {
    refTagCC_ += (**dit).name();
    if( i != 3 ) refTagCC_ += string(",");
  }
  // check if intermeidates or only four point diagrams
  bool intermediates(false);
  for(auto diagram : diagrams_) {
    if(diagram.intermediate) {
      intermediates=true;
      break;
    }
  }
  if(!intermediates) {
    incoming_= PDPtr();
    outgoing_.clear();
    generator()->log() << "Only four body diagrams for decay "
		       << refTag_  << " in GeneralThreeBodyDecayer::"
		       << "setDecayInfo(), omitting decay\n";
    return false;
  }
  // set the colour factors and return the answer
  if(setColourFactors(symfac)) return true;
  incoming_= PDPtr();
  outgoing_.clear();
  return false;
}

void GeneralThreeBodyDecayer::doinit() {
  DecayIntegrator::doinit();
  if(outgoing_.empty()) return;
  setupDiagrams(false);
}

void GeneralThreeBodyDecayer::doinitrun() {
  if(outgoing_.empty()) return;
  setupDiagrams(true);
  DecayIntegrator::doinitrun();
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
  if(!widthCalc_) {
    string tag = incoming_->name() + "->";
    tag += outgoing_[0]->name() + "," + outgoing_[1]->name() + ","
      + outgoing_[2]->name() + ";";
    DMPtr dm = generator()->findDecayMode(tag);
    widthCalc_ = threeBodyMEIntegrator(*dm);
  }
  return widthCalc_->partialWidth(sqr(inpart.second));
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
    switch(diagrams_[idiag].channelType) {
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
    else if(triplet.size()==3) {
      tColinePtr col[3] = {ColourLine::create(outgoing[0]),
			   ColourLine::create(outgoing[1]),
			   ColourLine::create(outgoing[2])};
      col[0]->setSourceNeighbours(col[1],col[2]);
    }
    else if(antitriplet.size()==3) {
      tColinePtr col[3] = {ColourLine::create(outgoing[0],true),
			   ColourLine::create(outgoing[1],true),
			   ColourLine::create(outgoing[2],true)};
      col[0]->setSinkNeighbours(col[1],col[2]);
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
    else if (singlet.size()==1&&antitriplet.size()==2) {
      tColinePtr col[2] = {ColourLine::create(outgoing[antitriplet[0]],true),
			   ColourLine::create(outgoing[antitriplet[1]],true)};
      parent.colourLine()->setSinkNeighbours(col[0],col[1]);
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
    else if (singlet.size()==1&&triplet.size()==2) {
      tColinePtr col[2] = {ColourLine::create(outgoing[triplet[0]]),
			   ColourLine::create(outgoing[triplet[1]])};
      parent.antiColourLine()->setSourceNeighbours(col[0],col[1]);
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
    else if(triplet.size()==3) {
      tColinePtr col[2];
      if(colourFlow()==0) {
	outgoing[0]->incomingColour    (const_ptr_cast<tPPtr>(&parent));
	col[0] = ColourLine::create(outgoing[1]);
	col[1] = ColourLine::create(outgoing[2]);
      }
      else if(colourFlow()==1) {
	outgoing[1]->incomingColour    (const_ptr_cast<tPPtr>(&parent));
	col[0] = ColourLine::create(outgoing[0]);
	col[1] = ColourLine::create(outgoing[2]);
      }
      else if(colourFlow()==2) {
	outgoing[2]->incomingColour    (const_ptr_cast<tPPtr>(&parent));
	col[0] = ColourLine::create(outgoing[0]);
	col[1] = ColourLine::create(outgoing[1]);
      }
      else
	assert(false);
      parent.antiColourLine()->setSourceNeighbours(col[0],col[1]);
    }
    else if(antitriplet.size()==3) {
      tColinePtr col[2];
      if(colourFlow()==0) {
	outgoing[0]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
	col[0] = ColourLine::create(outgoing[1],true);
	col[1] = ColourLine::create(outgoing[2],true);
      }
      else if(colourFlow()==1) {
	outgoing[1]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
	col[0] = ColourLine::create(outgoing[0],true);
	col[1] = ColourLine::create(outgoing[2],true);
      }
      else if(colourFlow()==2) {
	outgoing[2]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
	col[0] = ColourLine::create(outgoing[0],true);
	col[1] = ColourLine::create(outgoing[1],true);
      }
      else
	assert(false);
      parent.colourLine()->setSinkNeighbours(col[0],col[1]);
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
  // check if any intermediate photons
  bool hasPhoton=false;
  for(unsigned int iy=0;iy<diagmap_.size();++iy) {
    unsigned int ix=diagmap_[iy];
    if(getProcessInfo()[ix].intermediate->id()==ParticleID::gamma)
      hasPhoton = true;
  }
  // loop over channels
  Energy min = incoming()->mass();
  int nchannel(0);
  pair<int,Energy> imin[4]={make_pair(-1,-1.*GeV),make_pair(-1,-1.*GeV),
			    make_pair(-1,-1.*GeV),make_pair(-1,-1.*GeV)};
  Energy absmin = -1e20*GeV;
  int minType   = -1;
  for(unsigned int iy=0;iy<diagmap_.size();++iy) {
    unsigned int ix=diagmap_[iy];
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
    if((dm1<ZERO||dm2<ZERO)&&!hasPhoton) {
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
      inwidth.push_back(widthOption() ==3 ? ZERO : getProcessInfo()[ix].intermediate->width());
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
  else if(intOpt_==1&&minType>0&&getProcessInfo()[imin[minType].first].intermediate->id()!=ParticleID::gamma) {
    intype.push_back(minType);
    inpow.push_back(0.);
    inmass.push_back(getProcessInfo()[imin[minType].first].intermediate->mass());
    inwidth.push_back(widthOption() ==3 ? ZERO : getProcessInfo()[imin[minType].first].intermediate->width());
    inweights = vector<double>(1,1.);
    return;
  }
  for(unsigned int ix=1;ix<4;++ix) {
    if(imin[ix].first>=0) {
      intype.push_back(ix);
      if(getProcessInfo()[imin[ix].first].intermediate->id()!=ParticleID::gamma) {
	inpow.push_back(0.);
	inmass.push_back(getProcessInfo()[imin[ix].first].intermediate->mass());
	inwidth.push_back(widthOption() ==3 ? ZERO : getProcessInfo()[imin[ix].first].intermediate->width());
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

bool GeneralThreeBodyDecayer::setColourFactors(double symfac) {
  string name = incoming_->PDGName() + "->";
  vector<int> sng,trip,atrip,oct;
  unsigned int iloc(0);

  for(vector<PDPtr>::const_iterator it = outgoing_.begin();
      it != outgoing_.end();++it) {
    name += (**it).PDGName() + " ";
    if     ((**it).iColour() == PDT::Colour0    ) sng.push_back(iloc) ;
    else if((**it).iColour() == PDT::Colour3    ) trip.push_back(iloc) ;
    else if((**it).iColour() == PDT::Colour3bar ) atrip.push_back(iloc);
    else if((**it).iColour() == PDT::Colour8    ) oct.push_back(iloc) ;
    ++iloc;
  }
  // colour neutral decaying particle
  if     ( incoming_->iColour() == PDT::Colour0) {
    // options are all neutral or triplet/antitriplet+ neutral
    if(sng.size()==3) {
      nflow_ = 1;
      colour_        = vector<DVector>(1,DVector(1,1.));
      colourLargeNC_ = vector<DVector>(1,DVector(1,1.));
    }
    else if(sng.size()==1&&trip.size()==1&&atrip.size()==1) {
      nflow_ = 1;
      colour_         = vector<DVector>(1,DVector(1,3.));
      colourLargeNC_  = vector<DVector>(1,DVector(1,3.));
    }
    else if(trip.size()==1&&atrip.size()==1&&oct.size()==1) {
      nflow_ = 1;
      colour_         = vector<DVector>(1,DVector(1,4.));
      colourLargeNC_  = vector<DVector>(1,DVector(1,4.));
    }
    else if( trip.size() == 3 || atrip.size() == 3 ) {
      nflow_ = 1;
      colour_         = vector<DVector>(1,DVector(1,6.));
      colourLargeNC_  = vector<DVector>(1,DVector(1,6.));
      for(unsigned int ix=0;ix<diagrams_.size();++ix) {
	double sign = diagrams_[ix].channelType == TBDiagram::channel13 ? -1. : 1.;
	diagrams_[ix].       colourFlow = vector<CFPair>(1,make_pair(1,sign));
	diagrams_[ix].largeNcColourFlow = vector<CFPair>(1,make_pair(1,sign));
      }
    }
    else {
      generator()->log() << "Unknown colour flow structure for "
			 << "colour neutral decay " 
			 << name  << " in GeneralThreeBodyDecayer::"
			 << "setColourFactors(), omitting decay\n";
      return false;
    }
  }
  // colour triplet decaying particle
  else if( incoming_->iColour() == PDT::Colour3) {
    if(sng.size()==2&&trip.size()==1) {
      nflow_ = 1;
      colour_        = vector<DVector>(1,DVector(1,1.));
      colourLargeNC_ = vector<DVector>(1,DVector(1,1.));
    }
    else if(trip.size()==2&&atrip.size()==1) {
      nflow_ = 2;
      colour_.clear();
      colour_.resize(2,DVector(2,0.));
      colour_[0][0] = 3.; colour_[0][1] = 1.;
      colour_[1][0] = 1.; colour_[1][1] = 3.;
      colourLargeNC_.clear();
      colourLargeNC_.resize(2,DVector(2,0.));
      colourLargeNC_[0][0] = 3.; colourLargeNC_[1][1] = 3.;
      // sort out the contribution of the different diagrams to the colour
      // flows
      for(unsigned int ix=0;ix<diagrams_.size();++ix) {
	// colour singlet intermediate
	if(diagrams_[ix].intermediate->iColour()==PDT::Colour0) {
	  if(diagrams_[ix].channelType==trip[0]) {
	    diagrams_[ix].       colourFlow = vector<CFPair>(1,make_pair(1,1.));
	    diagrams_[ix].largeNcColourFlow = vector<CFPair>(1,make_pair(1,1.));
	  }
	  else {
	    diagrams_[ix].colourFlow        = vector<CFPair>(1,make_pair(2,1.));
	    diagrams_[ix].largeNcColourFlow = vector<CFPair>(1,make_pair(2,1.));
	  }
	}
	// colour octet intermediate
	else if(diagrams_[ix].intermediate->iColour()==PDT::Colour8) {
	  if(diagrams_[ix].channelType==trip[0]) {
	    vector<CFPair> flow(1,make_pair(2, 0.5  ));
	    diagrams_[ix].largeNcColourFlow = flow;
	    flow.push_back(       make_pair(1,-1./6.));
	    diagrams_[ix].colourFlow=flow;
	  }
	  else {
	    vector<CFPair> flow(1,make_pair(1, 0.5  ));
	    diagrams_[ix].largeNcColourFlow = flow;
	    flow.push_back(       make_pair(2,-1./6.));
	    diagrams_[ix].colourFlow=flow;
	  }
	}
	else {
	  generator()->log() << "Unknown colour for the intermediate in "
			     << "triplet -> triplet triplet antitriplet in "
			     << "GeneralThreeBodyDecayer::setColourFactors()"
			     << " for " << name << " omitting decay\n";
	  return false;
	}
      }
    }
    else if(trip.size()==1&&oct.size()==1&&sng.size()==1) {
      nflow_ = 1;
      colour_        = vector<DVector>(1,DVector(1,4./3.));
      colourLargeNC_ = vector<DVector>(1,DVector(1,4./3.));
    }
    else if(sng.size()==1&&atrip.size()==2) {
      nflow_ = 1;
      colour_        = vector<DVector>(1,DVector(1,2.));
      colourLargeNC_ = vector<DVector>(1,DVector(1,2.));
    }
    else {
      generator()->log() << "Unknown colour structure for "
			 << "triplet decay in "
			 << "GeneralThreeBodyDecayer::setColourFactors()"
			 << " for " << name << " omitting decay\n";
      return false;
    }
  }
  // colour antitriplet decaying particle
  else if( incoming_->iColour() == PDT::Colour3bar) {
    if(sng.size()==2&&atrip.size()==1) {
      nflow_ = 1;
      colour_        = vector<DVector>(1,DVector(1,1.));
      colourLargeNC_ = vector<DVector>(1,DVector(1,1.));
    }
    else if(atrip.size()==2&&trip.size()==1) {
      nflow_ = 2;
      colour_.clear();
      colour_.resize(2,DVector(2,0.));
      colour_[0][0] = 3.; colour_[0][1] = 1.;
      colour_[1][0] = 1.; colour_[1][1] = 3.;
      colourLargeNC_.clear();
      colourLargeNC_.resize(2,DVector(2,0.));
      colourLargeNC_[0][0] = 3.; colourLargeNC_[1][1] = 3.;
      // sort out the contribution of the different diagrams to the colour
      // flows
      for(unsigned int ix=0;ix<diagrams_.size();++ix) {
	// colour singlet intermediate
	if(diagrams_[ix].intermediate->iColour()==PDT::Colour0) {
	  if(diagrams_[ix].channelType==atrip[0]) {
	    diagrams_[ix].       colourFlow = vector<CFPair>(1,make_pair(1,1.));
	    diagrams_[ix].largeNcColourFlow = vector<CFPair>(1,make_pair(1,1.));
	  }
	  else {
	    diagrams_[ix].colourFlow        = vector<CFPair>(1,make_pair(2,1.));
	    diagrams_[ix].largeNcColourFlow = vector<CFPair>(1,make_pair(2,1.));
	  }
	}
	// colour octet intermediate
	else if(diagrams_[ix].intermediate->iColour()==PDT::Colour8) {
	  if(diagrams_[ix].channelType==atrip[0]) {
	    vector<CFPair> flow(1,make_pair(2, 0.5  ));
	    diagrams_[ix].largeNcColourFlow = flow;
	    flow.push_back(       make_pair(1,-1./6.));
	    diagrams_[ix].colourFlow=flow;
	  }
	  else {
	    vector<CFPair> flow(1,make_pair(1, 0.5  ));
	    diagrams_[ix].largeNcColourFlow = flow;
	    flow.push_back(       make_pair(2,-1./6.));
	    diagrams_[ix].colourFlow=flow;
	  }
	}
	else {
	  generator()->log() << "Unknown colour for the intermediate in "
			     << "antitriplet -> antitriplet antitriplet triplet in "
			     << "GeneralThreeBodyDecayer::setColourFactors()"
			     << " for " << name << " omitting decay\n";
	  return false;
	}
      }
    }
    else if(atrip.size()==1&&oct.size()==1&&sng.size()==1) {
      nflow_ = 1;
      colour_        = vector<DVector>(1,DVector(1,4./3.));
      colourLargeNC_ = vector<DVector>(1,DVector(1,4./3.));
    }
    else if(sng.size()==1&&trip.size()==2) {
      nflow_ = 1;
      colour_        = vector<DVector>(1,DVector(1,2.));
      colourLargeNC_ = vector<DVector>(1,DVector(1,2.));
    }
    else {
      generator()->log() << "Unknown colour antitriplet decay in "
			 << "GeneralThreeBodyDecayer::setColourFactors()"
			 << " for " << name << " omitting decay\n";
      return false;
    }
  }
  // colour octet particle
  else if( incoming_->iColour() == PDT::Colour8) {
    // triplet antitriplet
    if(trip.size() == 1 && atrip.size() == 1 && sng.size() == 1) {
      nflow_ = 1;
      colour_        = vector<DVector>(1,DVector(1,0.5));
      colourLargeNC_ = vector<DVector>(1,DVector(1,0.5));
    }
    // three (anti)triplets
    else if(trip.size()==3||atrip.size()==3) {
      nflow_ = 3;
      colour_        = vector<DVector>(3,DVector(3,0.));
      colourLargeNC_ = vector<DVector>(3,DVector(3,0.));
      colour_[0][0] = 1.; colour_[1][1] = 1.; colour_[2][2] = 1.;
      colour_[0][1] = -0.5; colour_[1][0] = -0.5;
      colour_[0][2] = -0.5; colour_[2][0] = -0.5;
      colour_[1][2] = -0.5; colour_[2][1] = -0.5;
      colourLargeNC_ = vector<DVector>(3,DVector(3,0.));
      colourLargeNC_[0][0] = 1.; colourLargeNC_[1][1] = 1.; colourLargeNC_[2][2] = 1.;
      // sett the factors for the diagrams
      for(unsigned int ix=0;ix<diagrams_.size();++ix) {
	tPDPtr inter = diagrams_[ix].intermediate;
	if(inter->CC()) inter = inter->CC();
	unsigned int io[2]={1,2};
	double sign = diagrams_[ix].channelType == TBDiagram::channel13 ? -1. : 1.;
	for(unsigned int iy=0;iy<3;++iy) {
	  if     (iy==1) io[0]=0;
	  else if(iy==2) io[1]=1;
	  tPDVector decaylist = diagrams_[ix].vertices.second->search(iy, inter);
	  if(decaylist.empty()) continue;
	  bool found=false;
	  for(unsigned int iz=0;iz<decaylist.size();iz+=3) {	    
	    if(decaylist[iz+io[0]]->id()==diagrams_[ix].outgoingPair.first &&
	       decaylist[iz+io[1]]->id()==diagrams_[ix].outgoingPair.second) {
	      sign *= 1.;
	      found = true;
	    }
	    else if(decaylist[iz+io[0]]->id()==diagrams_[ix].outgoingPair.second &&
		    decaylist[iz+io[1]]->id()==diagrams_[ix].outgoingPair.first ) {
	      sign *= -1.;
	      found = true;
	    }
	  }
	  if(found) {
	    if(iy==1) sign *=-1.;
	    break;
	  }
	}
	diagrams_[ix].       colourFlow = vector<CFPair>(1,make_pair(diagrams_[ix].channelType+1,sign));
	diagrams_[ix].largeNcColourFlow = vector<CFPair>(1,make_pair(diagrams_[ix].channelType+1,sign));
      }
    }
    // unknown
    else {
      generator()->log() << "Unknown colour octet decay in "
			 << "GeneralThreeBodyDecayer::setColourFactors()"
			 << " for " << name << " omitting decay\n";
      return false;
    }
  }
  else if (incoming_->iColour() == PDT::Colour6 ) {
    generator()->log() << "Unknown colour sextet decay in "
		       << "GeneralThreeBodyDecayer::setColourFactors()"
		       << " for " << name << " omitting decay\n";
    return false;
  }
  else if (incoming_->iColour() == PDT::Colour6bar ) {
    generator()->log() << "Unknown colour anti-sextet decay in "
		       << "GeneralThreeBodyDecayer::setColourFactors()"
		       << " for " << name << " omitting decay\n";
    return false;
  }

  assert(nflow_ != 999);

  for(unsigned int ix=0;ix<nflow_;++ix) {
    for(unsigned int iy=0;iy<nflow_;++iy) {
      colour_       [ix][iy] /= symfac;
      colourLargeNC_[ix][iy] /= symfac;
    }
  }
  if( Debug::level > 1 ) {
    generator()->log() << "Mode: " << name << " has colour factors\n";
    for(unsigned int ix=0;ix<nflow_;++ix) {
      for(unsigned int iy=0;iy<nflow_;++iy) {
	generator()->log() << colour_[ix][iy] << " ";
      }
      generator()->log() << "\n";
    }
    for(unsigned int ix=0;ix<diagrams_.size();++ix) {
      generator()->log() << "colour flow for diagram : " << ix;
      for(unsigned int iy=0;iy<diagrams_[ix].colourFlow.size();++iy)
	generator()->log() << "(" << diagrams_[ix].colourFlow[iy].first  << "," 
			   << diagrams_[ix].colourFlow[iy].second << "); ";
      generator()->log() << "\n";
    }
  }
  return true;
}

void GeneralThreeBodyDecayer::setupDiagrams(bool kinCheck) {
  clearModes();
  // create the phase space integrator
  tPDVector extpart(1,incoming_);
  extpart.insert(extpart.end(),outgoing_.begin(),outgoing_.end());
  // create the integration channels for the decay
  DecayPhaseSpaceModePtr mode(new_ptr(DecayPhaseSpaceMode(extpart,this,true)));
  DecayPhaseSpaceChannelPtr newchannel;
  // create the phase-space channels for the integration
  unsigned int nmode(0);
  unsigned int idiag(0);
  for(vector<TBDiagram>::iterator it = diagrams_.begin();it!=diagrams_.end();++it) {
    if(it->channelType==TBDiagram::fourPoint||
       it->channelType==TBDiagram::UNDEFINED) {
      idiag+=1;
      continue;
    }
    // create the new channel
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    int jac = 0;
    double power = 0.0;
    if ( it->intermediate->mass() == ZERO ||
  	 it->intermediate->width() == ZERO ) {
      jac = 1;
      power = -2.0;
    }
    if(it->channelType==TBDiagram::channel23) {
      newchannel->addIntermediate(extpart[0],0,0.0,-1,1);
      newchannel->addIntermediate(it->intermediate,jac,power, 2,3);
    }
    else if(it->channelType==TBDiagram::channel13) {
      newchannel->addIntermediate(extpart[0],0,0.0,-1,2);
      newchannel->addIntermediate(it->intermediate,jac,power, 1,3);
    }
    else if(it->channelType==TBDiagram::channel12) {
      newchannel->addIntermediate(extpart[0],0,0.0,-1,3);
      newchannel->addIntermediate(it->intermediate,jac,power, 1,2);
    }
    newchannel->init();
    if(kinCheck&&!newchannel->checkKinematics()) {
      generator()->log() << "Erasing diagram "
			 << *it
			 << "from three body decay as zero width propagator can be on-shell,\n"
			 << "hopefully this diagram is zero in your model, but you should check this\n";
      it = diagrams_.erase(it);
      if(it == diagrams_.end()) break;
      continue;
    }
    diagmap_.push_back(idiag);
    mode->addChannel(newchannel);
    ++nmode;
    ++idiag;
  }
  if(nmode==0) {
    string mode = extpart[0]->PDGName() + "->";
    for(unsigned int ix=1;ix<extpart.size();++ix) mode += extpart[ix]->PDGName() + " ";
    generator()->log() << "No decay channels in GeneralThreeBodyDecayer::"
		       << "setupDiagrams() for " << mode << "\n";
    return;
  }
  // add the mode
  vector<double> wgt(nmode,1./double(nmode));
  addMode(mode,1.,wgt);
}
