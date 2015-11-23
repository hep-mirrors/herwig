// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralFourBodyDecayer class.
//

#include "GeneralFourBodyDecayer.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
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
     "of all four-body decays based on spin structures in Herwig.");

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
    if( i != 4 ) testtag += string(",");
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
					  const vector<NBDiagram> & process,
					  double symfac) {
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
    if( dit != refmode.begin() )  _reftagcc += string(",");
    _reftagcc += (**dit).name();
  }
  // set the colour factors and return the answer
  return setColourFactors(symfac);
}

Energy GeneralFourBodyDecayer::partialWidth(tPDPtr inpart,
					    OrderedParticles outgoing) const {
  tPDVector temp = tPDVector(outgoing.begin(),outgoing.end());
  bool cc=false;
  int imode = modeNumber(cc,inpart,temp);
  if(imode<0) return ZERO;
  else return initializePhaseSpaceMode(0,true,true);  
}

namespace {
  bool addIntermediate(DecayPhaseSpaceChannelPtr newChannel,
		       const NBVertex &vertex,vector<unsigned int> &order,
		       int & ioff, int & loc) {
    if(vertex.vertices.size()!=2) return false;
    assert(!vertex.vertices.begin()->second.incoming);
    int first = order[loc];
    ++loc;
    unsigned int type = vertex.incoming->width()!=ZERO ? 0 : 1;
    double      power = vertex.incoming->width()!=ZERO ? 0. : -2.;
    // one off-shell
    if((++vertex.vertices.begin())->second.incoming) {
      newChannel->addIntermediate(vertex.incoming,type,power,ioff,first);
      --ioff;
      return addIntermediate(newChannel,(++vertex.vertices.begin())->second,
			     order,ioff,loc);
    }
    // no off-shell
    else {
      int second = order[loc];
      ++loc;
      newChannel->addIntermediate(vertex.incoming,type,power,first,second);
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
    int ioff = -1;
    int loc=0;
    // two off-shell particles
    if(_diagrams[ix].vertices.begin()->second.incoming) {
      newChannel->addIntermediate(extpart[0],0,0.0,-1,-2);
      if(!addIntermediate(newChannel,(  _diagrams[ix].vertices.begin())->second,
			  _diagrams[ix].channelType,ioff,loc)) continue;
      if(!addIntermediate(newChannel,(++_diagrams[ix].vertices.begin())->second,
			  _diagrams[ix].channelType,ioff,loc)) continue;
    }
    // only one off-shell particle
    else {
      int first = _diagrams[ix].channelType[loc];
      ++loc;
      newChannel->addIntermediate(extpart[0],0,0.0,-1,first);
      --ioff;
      if(!addIntermediate(newChannel,(++_diagrams[ix].vertices.begin())->second,
			  _diagrams[ix].channelType,ioff,loc)) continue;
    }
    _diagmap.push_back(ix);
    mode->addChannel(newChannel);
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
    assert(false);
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
      throw Exception() << "Unknown colour for particle in GeneralFourBodyDecayer::"
			<< "colourConnections()" << Exception::runerror;
    }
  }
  // colour neutral decaying particle
  if     ( inColour == PDT::Colour0) {
    // options are all neutral
    if(singlet.size()==4) return;
    // two singlets + triplet/antitriplet
    else if(singlet.size()==2&&triplet.size()==1&&antitriplet.size()==1) {
      outgoing[triplet[0]]->antiColourNeighbour(outgoing[antitriplet[0]]);
    }
    // 2 x triplet + antitriplet
    else if(triplet.size()==2&&antitriplet.size()==2) {
      if(colourFlow()==0) {
	outgoing[triplet[0]]->antiColourNeighbour(outgoing[antitriplet[0]]);
	outgoing[triplet[1]]->antiColourNeighbour(outgoing[antitriplet[1]]);
      }
      else {
	outgoing[triplet[0]]->antiColourNeighbour(outgoing[antitriplet[1]]);
	outgoing[triplet[1]]->antiColourNeighbour(outgoing[antitriplet[0]]);
      }
    }
    // 3 triplets
    else if(triplet.size()==3&&singlet.size()==1&&antitriplet.size()==0) {
      tColinePtr col[3] = {ColourLine::create(outgoing[triplet[0]]),
			   ColourLine::create(outgoing[triplet[1]]),
			   ColourLine::create(outgoing[triplet[2]])};
      col[0]->setSourceNeighbours(col[1],col[2]);
    }
    // 3 antitriplets
    else if(triplet.size()==0&&singlet.size()==1&&antitriplet.size()==3) {
      tColinePtr col[3] = {ColourLine::create(outgoing[antitriplet[0]],true),
			   ColourLine::create(outgoing[antitriplet[1]],true),
			   ColourLine::create(outgoing[antitriplet[2]],true)};
      col[0]->setSinkNeighbours(col[1],col[2]);
    }
    else {
      string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
	+ out[1]->PDGName() + " " + out[2]->PDGName()+ " " + out[3]->PDGName();
      throw Exception() 
	<< "Unknown colour structure in GeneralFourBodyDecayer::"
	<< "colourConnections() for singlet decaying particle "
	<< mode << Exception::runerror;
    }
  }
  // colour triplet decaying particle
  else if( inColour == PDT::Colour3) {
    if(singlet.size()==3&&triplet.size()==1) {
      outgoing[triplet[0]]->incomingColour(const_ptr_cast<tPPtr>(&parent));
    }
    else if(antitriplet.size()==1&&triplet.size()==2&&singlet.size()==1) {
      if(colourFlow()==0) {
	outgoing[triplet[0]]->incomingColour(const_ptr_cast<tPPtr>(&parent));
	outgoing[antitriplet[0]]->colourNeighbour(outgoing[triplet[1]]);
      }
      else {
	outgoing[triplet[1]]->incomingColour(const_ptr_cast<tPPtr>(&parent));
	outgoing[antitriplet[0]]->colourNeighbour(outgoing[triplet[0]]);
      }
    }
    else {
      string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
	+ out[1]->PDGName() + " " + out[2]->PDGName() + " " + out[3]->PDGName();
      throw Exception() 
	<< "Unknown colour structure in GeneralFourBodyDecayer::"
	<< "colourConnections() for triplet decaying particle " 
	<< mode << Exception::runerror;
    }
  }
  else if( inColour == PDT::Colour3bar) {
    if(singlet.size()==2&&antitriplet.size()==1) {
      outgoing[antitriplet[0]]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
    }
    else if(antitriplet.size()==2&&triplet.size()==1&&singlet.size()==1) {
      if(colourFlow()==0) {
	outgoing[antitriplet[0]]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
	outgoing[triplet[0]]->antiColourNeighbour(outgoing[antitriplet[1]]);
      }
      else {
	outgoing[antitriplet[1]]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
	outgoing[triplet[0]]->antiColourNeighbour(outgoing[antitriplet[0]]);
      }
    }
    else {
      string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
	+ out[1]->PDGName() + " " + out[2]->PDGName() + " " + out[3]->PDGName();
      throw Exception() 
	<< "Unknown colour structure in GeneralFourBodyDecayer::"
	<< "colourConnections() for anti-triplet decaying particle" 
	<< mode << Exception::runerror;
    }
  }
  else {
    string mode = parent.PDGName() + " -> " + out[0]->PDGName() + " "
      + out[1]->PDGName() + " " + out[2]->PDGName()+ " " + out[3]->PDGName();
    throw Exception() 
      << "Unknown colour for decaying particle in GeneralFourBodyDecayer::"
      << "colourConnections() "
      << mode << Exception::runerror;
  }
}

bool GeneralFourBodyDecayer::setColourFactors(double symfac) {
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
    else if(trip.size()==2&&atrip.size()==2) {
      _nflow = 2;
      _colour.clear();
      _colour.resize(2,DVector(2,0.));
      _colour[0][0] = 9.; _colour[0][1] = 3.;
      _colour[1][0] = 3.; _colour[1][1] = 9.;
      _colourLargeNC.clear();
      _colourLargeNC.resize(2,DVector(2,0.));
      _colourLargeNC[0][0] = 9.; _colourLargeNC[1][1] = 9.;
      // find potential colour connected pairs
      pair<unsigned int,unsigned int> first(make_pair(0,0)),second(make_pair(0,0));
      for(unsigned int ix=0;ix<_outgoing.size();++ix) {
 	if(_outgoing[ix]->iColour()==PDT::Colour3) {
	  if(first.first ==0) first .first  = ix+1;
	  else                second.first  = ix+1;
	}
	else if(_outgoing[ix]->iColour()==PDT::Colour3bar) {
	  if(first.second==0) first .second = ix+1;
	  else                second.second = ix+1;
	}
      }
      // sort out the contribution of the different diagrams
      // to the colour flows
      for(vector<NBDiagram>::iterator it = _diagrams.begin();
	  it!=_diagrams.end();++it) {
	unsigned int iloc[4]={ it->channelType[0], it->channelType[1],
			       it->channelType[2], it->channelType[3]};
	if(_outgoing[iloc[1]-1]->iColour()==PDT::Colour3) swap(iloc[0],iloc[1]);
	if(_outgoing[iloc[3]-1]->iColour()==PDT::Colour3) swap(iloc[2],iloc[3]);
	if(iloc[0]>iloc[2]) {
	  swap(iloc[0],iloc[2]);
	  swap(iloc[1],iloc[3]);
	}
	assert(iloc[0]==first .first);
	assert(iloc[2]==second.first);
	tPDPair inter;
	if(it->vertices.begin()->second.incoming) {
	  inter.first  =     it->vertices.begin()  ->second.incoming;
	  inter.second = (++(it->vertices.begin()))->second.incoming;
	}
	else {
	  const NBVertex & vertex =  (++(it->vertices.begin()))->second;
	  inter.first  =  (++(it->vertices.begin()))->second.incoming;
	  inter.second = (++vertex.vertices.begin())->second.incoming;
	}
	// first topo two neutral, or second triplet + neutral
	if((inter. first->iColour()==PDT::Colour0 &&
	    inter.second->iColour()==PDT::Colour0) ||
	   (inter. first->iColour()==PDT::Colour3 &&
	    inter.second->iColour()==PDT::Colour0) ||
	   (inter. first->iColour()==PDT::Colour3bar &&
	    inter.second->iColour()==PDT::Colour0)) {
	  if(iloc[0]==first.second) {
	    it->       colourFlow = vector<CFPair>(1,make_pair(1,1.));
	    it->largeNcColourFlow = vector<CFPair>(1,make_pair(1,1.));
	  }
	  else {
	    it->colourFlow        = vector<CFPair>(1,make_pair(2,1.));
	    it->largeNcColourFlow = vector<CFPair>(1,make_pair(2,1.));
	  }
	}
	// first typo two octet, or second triplet + octet
	else if ((inter. first->iColour()==PDT::Colour8 &&
		  inter.second->iColour()==PDT::Colour8) ||
		 (inter. first->iColour()==PDT::Colour3 &&
		  inter.second->iColour()==PDT::Colour8) ||
		 (inter. first->iColour()==PDT::Colour3bar &&
		  inter.second->iColour()==PDT::Colour8)) {
	  if(iloc[0]==first.second) {
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
	  assert(false);
      }
    }
    else if ( trip.size() == 3 || atrip.size() == 3 ) {
      _nflow = 1;
      _colour         = vector<DVector>(1,DVector(1,6.));
      _colourLargeNC  = vector<DVector>(1,DVector(1,6.));
      vector<int> coloured = trip.size()==3 ? trip : atrip;
      for(vector<NBDiagram>::iterator it = _diagrams.begin();
	  it!=_diagrams.end();++it) {
	vector<int> coloured = trip.size()==3 ? trip : atrip;
	// get the ordering of the triplets
	int iloc[4];
	for(unsigned int ix=0;ix<4;++ix) iloc[ix] = it->channelType[ix]-1;
	double sign(1.);
	if(iloc[coloured[0]]>iloc[coloured[2]]) {
	  swap(coloured[0],coloured[2]);
	  sign *= -1.;
	}
	if(iloc[coloured[0]]>iloc[coloured[1]]) {
	  swap(coloured[0],coloured[1]);
	  sign *= -1.;
	}
	if(iloc[coloured[1]]>iloc[coloured[2]]) {
	  swap(coloured[1],coloured[2]);
	  sign *= -1.;
	}
	// extract the vertices
	// get the other vertices
	const NBVertex & second = (*it).vertices.begin()->second.incoming ?
	  (*it).vertices.begin()->second : (++(*it).vertices.begin())->second;
	// get the other vertices
	const NBVertex & third = (*it).vertices.begin()->second.incoming ?
	  (++(*it).vertices.begin())->second : (++second.vertices.begin())->second;
	// find the BV vertex
	VertexBasePtr vertex;
	tcPDPtr part[3];
	// first topologyw
	if(it->vertices.begin()->second.incoming) {
	  if(_outgoing[iloc[0]]->coloured()&&
	     _outgoing[iloc[1]]->coloured()) {
	    part[0] = _outgoing[iloc[0]];
	    part[1] = _outgoing[iloc[1]];
	    part[2] = second.incoming->CC() ? second.incoming->CC() : second.incoming;
	    vertex = second.vertex;
	  }
	  else {
	    part[0] = third .incoming->CC() ? third .incoming->CC() : third .incoming;
	    part[1] = _outgoing[iloc[2]];
	    part[2] = _outgoing[iloc[3]];
	    vertex = third.vertex;
	  }
	}
	else {
	  if(_outgoing[iloc[2]]->coloured()&&
	     _outgoing[iloc[3]]->coloured()) {
	    part[0] = third .incoming->CC() ? third .incoming->CC() : third .incoming;
	    part[1] = _outgoing[iloc[2]];
	    part[2] = _outgoing[iloc[3]];
	    vertex = third.vertex;
	  }
	  else {
	    part[0] = second.incoming->CC() ? second.incoming->CC() : second.incoming;
	    part[1] = _outgoing[iloc[1]];
	    part[2] = third .incoming->CC() ? third .incoming->CC() : third .incoming;
	    vertex = second.vertex;
	  }
	}
      	unsigned int io[2]={1,2};
      	for(unsigned int iy=0;iy<3;++iy) {
      	  if     (iy==1) io[0]=0;
      	  else if(iy==2) io[1]=1;
	  tPDVector decaylist = vertex->search(iy,part[0]);
      	  if(decaylist.empty()) continue;
      	  bool found=false;
       	  for(unsigned int iz=0;iz<decaylist.size();iz+=3) {	    
      	    if(decaylist[iz+io[0]]==part[1]&&
      	       decaylist[iz+io[1]]==part[2]) {
      	      sign *= 1.;
      	      found = true;
      	    }
      	    else if(decaylist[iz+io[0]]==part[2] &&
      		    decaylist[iz+io[1]]==part[1] ) {
      	      sign *= -1.;
      	      found = true;
      	    }
      	  }
      	  if(found) {
      	    if(iy==1) sign *=-1.;
      	    break;
      	  }
	}
	it->       colourFlow = vector<CFPair>(1,make_pair(1,sign));
	it->largeNcColourFlow = vector<CFPair>(1,make_pair(1,sign));
      }
    }
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
      _colour.clear();
      _colour.resize(2,DVector(2,0.));
      _colour[0][0] = 3.; _colour[0][1] = 1.;
      _colour[1][0] = 1.; _colour[1][1] = 3.;
      _colourLargeNC.clear();
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
	    if( it->channelType[0] == ifirst ||
		it->channelType[1] == ifirst) {
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
	    if( it->channelType[2] == ifirst ||
	        it->channelType[3] == ifirst) {
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
	    if( it->channelType[0] == ifirst ||
		it->channelType[1] == ifirst) {
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
	    if( it->channelType[2] == ifirst ||
	        it->channelType[3] == ifirst) {
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
	    iloc = it->channelType[1];
	    const NBVertex & vertex =  (++(it->vertices.begin()))->second;
	    inter = (++vertex.vertices.begin())->second.incoming;
	  }
	  else {
	    iloc = it->channelType[0];
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
      _colour.clear();
      _colour.resize(2,DVector(2,0.));
      _colour[0][0] = 3.; _colour[0][1] = 1.;
      _colour[1][0] = 1.; _colour[1][1] = 3.;
      _colourLargeNC.clear();
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
	    if( it->channelType[0] == ifirst ||
		it->channelType[1] == ifirst) {
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
	    if( it->channelType[2] == ifirst ||
	        it->channelType[3] == ifirst) {
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
	    if( it->channelType[0]       == ifirst ||
	       it->channelType[1]  == ifirst) {
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
	    if(it->channelType[2] == ifirst ||
	        it->channelType[3]      == ifirst) {
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
	    iloc = it->channelType[1];
	    const NBVertex & vertex =  (++(it->vertices.begin()))->second;
	    inter = (++vertex.vertices.begin())->second.incoming;
	  }
	  else {
	    iloc = it->channelType[0];
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
  for(unsigned int ix=0;ix<_nflow;++ix) {
    for(unsigned int iy=0;iy<_nflow;++iy) {
      _colour       [ix][iy] /= symfac;
      _colourLargeNC[ix][iy] /= symfac;
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




















