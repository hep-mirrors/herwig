// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HwRemDecayer class.
//

#include "HwRemDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <ThePEG/Interface/Reference.h>  
#include <ThePEG/Interface/Switch.h>  
#include "ThePEG/PDT/StandardMatchers.h"
#include "Herwig++/Shower/ShowerHandler.h"

using namespace Herwig;
using namespace ThePEG;

HwRemDecayer::~HwRemDecayer() {}

int HwRemDecayer::HadronContent::getValence() {
  if(extracted != -1)
    throw Exception() << "Try to extract second valence quark in "
		      << "HwRemDecayer::GetValence()"
		      << Exception::runerror; 
  int index( UseRandom::irnd( flav.size() ) );
  extracted = index;
  return sign*flav[index];
}

void HwRemDecayer::HadronContent::extract(int id) {
  for(unsigned int i=0; i<flav.size(); i++)
    if(id == sign*flav[i]){
      extracted = i;
      break;
    }      
}

bool HwRemDecayer::HadronContent::isSeaQuark(tcPPtr parton) const{
  return ((parton->id() != ParticleID::g) && 
	  ( !isValenceQuark(parton) ) );
}

bool HwRemDecayer::HadronContent::isValenceQuark(tcPPtr parton) const{
  int id(parton->id());
  for(unsigned int i=0; i<flav.size(); i++)
    if(id == sign*flav[i])
      return true;

  return false;
}

long HwRemDecayer::HadronContent::RemID() const{
  if(extracted == -1)
    throw Exception() << "Try to build a Diquark id without "
		      << "having extracted something in "
		      << "HwRemDecayer::RemID(...)"
		      << Exception::runerror; 
  if(flav.size()==2)//the hadron was a meson
    return sign*flav[(extracted+1)%2];
 
  long remId;
  int id1(sign*flav[(extracted+1)%3]), 
    id2(sign*flav[(extracted+2)%3]),
    sign(0), spin(0);

  if (abs(id1) > abs(id2)) swap(id1, id2);                                                                   
  sign = (id1 < 0) ? -1 : 1; // Needed for the spin 0/1 part                                                 
  remId = id2*1000+id1*100;                                                                                  
                                                                                                           
  // Now decide if we have spin 0 diquark or spin 1 diquark                                                  
  if(id1 == id2) spin = 3; // spin 1                                                                         
  else spin = 1; // otherwise spin 0                                                                         

  remId += sign*spin;  
  return remId;
}

HwRemDecayer::HadronContent
HwRemDecayer::getHadronContent(tcPPtr hadron) const {
  long id(hadron->id());
  if(abs(id) < 99)
    throw Exception() << "No hadron as beam particle "
		      << "HwRemDecayer::GetHadronContent(...)"
		      << Exception::runerror; 
  HadronContent hc;
  hc.sign = id < 0? -1: 1;
  hc.flav.push_back((id = abs(id)/10)%10);
  hc.flav.push_back((id /= 10)%10);
  hc.flav.push_back((id /= 10)%10);
  hc.extracted = -1;
  return hc;
}

bool HwRemDecayer::accept(const DecayMode &) const {
  return true;
}

bool HwRemDecayer::multiCapable() const {  
  return true;
} 

bool HwRemDecayer::
canHandle(tcPDPtr particle, tcPDPtr parton) const {
  return ( BaryonMatcher::Check(*particle) || MesonMatcher::Check(*particle) ) &&
    StandardQCDPartonMatcher::Check(*parton);
}

void HwRemDecayer::initialize(pair<tRemPPtr, tRemPPtr> rems, Step & step) {
  if(!theSplittingOnOff) return;
  tcPPair beam(generator()->currentEventHandler()->currentCollision()->incoming());

  thestep = &step;
  theContent.first = getHadronContent(beam.first);
  theContent.second = getHadronContent(beam.second);
  theUsed.first = Lorentz5Momentum();
  theUsed.second = Lorentz5Momentum();
  theMaps.first.clear();
  theMaps.second.clear();
  theX.first = 0.0;
  theX.second = 0.0;
  theRems = rems;

  if( parent(theRems.first) != beam.first ||
      parent(theRems.second) != beam.second )
    throw Exception() << "Remnant order wrong in "
		      << "HwRemDecayer::initialize(...)"
		      << Exception::runerror;

  return;
}

void HwRemDecayer::split(tPPtr parton, HadronContent & content, 
			 tRemPPtr rem, Lorentz5Momentum & used, 
			 PartnerMap &partners, bool first) {
  tcPPtr beam(parent(rem));

  if(rem==theRems.first)
    theX.first += parton->momentum().rho()/beam->momentum().rho();
  else
    theX.second += parton->momentum().rho()/beam->momentum().rho();

  double check;
  if(rem==theRems.first) check = theX.first;
  else check = theX.second;

  if(1.0-check < 1e-3)     
    throw ShowerHandler::ExtraScatterVeto();

  bool anti;
  Lorentz5Momentum lastp(parton->momentum());
  int lastID(parton->id());
  PPtr newSea, newValence;
  ColinePtr cl;
  //set the beam object to access the PDF's.
  theForcedSplitter->setBeam(beam);

  Energy oldQ(theForcedSplitter->getQspac());
  
  //do nothing if already valence quark
  if(first && content.isValenceQuark(parton)){ 
    //set the extracted value, because otherwise no RemID could be generated.
    content.extract(lastID);
    //add the particle to the colour partners
    partners.push_back(make_pair(parton, tPPtr()));
    //set the sign
    anti = parton->hasAntiColour();
    if(rem==theRems.first)
      theanti.first = anti;
    else
      theanti.second = anti;
    return;
  }
  //or gluon for secondaries
  if(!first && lastID == ParticleID::g){
    partners.push_back(make_pair(parton, tPPtr()));
    return; 
  }
  if( lastID != ParticleID::g ){
    //do the gluon splitting
    // Create the new parton with its momentum and parent/child relationship set
    if(rem==theRems.first)
      newSea = theForcedSplitter->
	forceSplit(rem, -lastID, oldQ, theX.first, 
		   lastp, used,1, thestep, first);
    else
      newSea = theForcedSplitter->
	forceSplit(rem, -lastID, oldQ, theX.second, 
		   lastp, used,1, thestep, first);

    cl = new_ptr(ColourLine());
    if(newSea->id() > 0) cl->addColoured(newSea);
    else cl->addAntiColoured(newSea);

  }

  if(!first){
    partners.push_back(make_pair(parton, newSea));
    return;
  }
  if( !content.isValenceQuark(parton) ){
    //final valence splitting

    //This was in the old code, but probably accidental:
    //oldQ = theForcedSplitter->getQspac();
    if(rem==theRems.first)
      newValence = theForcedSplitter->
	forceSplit(rem, content.getValence(), 
		   oldQ, theX.first, lastp, used, 2, thestep, first);
    else
      newValence = theForcedSplitter->
	forceSplit(rem, content.getValence(), 
		   oldQ, theX.second, lastp, used, 2, thestep, first);

    //case of a gluon going into the hard subprocess
    if( lastID == ParticleID::g ){
      partners.push_back(make_pair(parton, tPPtr()));

      anti = newValence->hasAntiColour();
      if(rem==theRems.first)
	theanti.first = anti;
      else
	theanti.second = anti;

      parton->colourLine(!anti)->addColoured(newValence, anti);
      return;
    }
    //The valence quark will always be connected to the sea quark with opposite sign
    tcPPtr particle;
    if(lastID*newValence->id() < 0){
      particle = parton;
      partners.push_back(make_pair(newSea, tPPtr()));
    }else{
      particle = newSea;
      partners.push_back(make_pair(parton, tPPtr()));
    }
    anti = newValence->hasAntiColour();
    if(rem==theRems.first)
      theanti.first = anti;
    else
      theanti.second = anti;

    if(particle->colourLine()) particle->colourLine()->addAntiColoured(newValence);
    if(particle->antiColourLine()) particle->antiColourLine()->addColoured(newValence);  
    
  }
  return;
}

void HwRemDecayer::setRemMasses() const {
  // get the masses of the remnants
  Energy mrem[2];
  Lorentz5Momentum ptotal,pnew[2];
  vector<tRemPPtr> theprocessed;
  theprocessed.push_back(theRems.first);
  theprocessed.push_back(theRems.second);

  for(unsigned int ix=0;ix<2;++ix) {
    pnew[ix]=Lorentz5Momentum();
    for(unsigned int iy=0;iy<theprocessed[ix]->children().size();++iy) {
      pnew[ix]+=theprocessed[ix]->children()[iy]->momentum();
    }
    mrem[ix]=sqrt(pnew[ix].m2());
  }
  // now find the remnant remnant cmf frame
  Lorentz5Momentum prem[2]={theprocessed[0]->momentum(),
			    theprocessed[1]->momentum()};
  ptotal=prem[0]+prem[1];
  ptotal.rescaleMass();

  // boost momenta to this frame
  if(ptotal.m()< (pnew[0].m()+pnew[1].m())) throw Exception() 
    << "Not enough energy in both remnants in " 
    << "HwRemDecayer::setRemMasses() " 
    << Exception::eventerror;

  Boost boostv(-ptotal.boostVector());
  ptotal.boost(boostv);
  for(unsigned int ix=0;ix<2;++ix) {
    prem[ix].boost(boostv);
    // set the masses and energies,
    prem[ix].setMass(mrem[ix]);
    prem[ix].setE(0.5/ptotal.m()*(sqr(ptotal.m())+sqr(mrem[ix])-sqr(mrem[1-ix])));
    prem[ix].rescaleRho();
    // boost back to the lab
    prem[ix].boost(-boostv);
    // set the momenta of the remnants
    theprocessed[ix]->set5Momentum(prem[ix]);
  }

  // boost the decay products
  Lorentz5Momentum ptemp; 
  for(unsigned int ix=0;ix<2;++ix) { 
    Boost btorest(-pnew[ix].boostVector()); 
    Boost bfmrest( prem[ix].boostVector()); 
    for(unsigned int iy=0;iy<theprocessed[ix]->children().size();++iy) { 
      ptemp=theprocessed[ix]->children()[iy]->momentum(); 
      ptemp.boost(btorest); 
      ptemp.boost(bfmrest); 
      theprocessed[ix]->children()[iy]->set5Momentum(ptemp); 
    } 
  }
}

void HwRemDecayer::fixColours(PartnerMap partners, bool anti) const {
  PartnerMap::const_iterator prev;
  tPPtr pnew, pold;
  ColinePtr clnew, clold;

  assert(partners.size()>=2);
  for(PartnerMap::iterator it=partners.begin(); 
      it!=partners.end(); it++){
    if(it==partners.begin()) continue;
    prev = it - 1;

    //determine the particles to work with
    pold = prev->first;
    if(prev->second){
      if(pold->hasAntiColour() != anti)
	pold = prev->second;
    }
    assert(pold);

    pnew = it->first;
    if(it->second){
      if(it->second->colourLine(!anti)) //look for the opposite colour
	pnew = it->second;
    }
    assert(pnew);


    //save the corresponding colour lines
    clold = pold->colourLine(anti);
    clnew = pnew->colourLine(!anti);

    assert(clold);

    
    if(clnew){//there is already a colour line (not the final diquark)

      if( (clnew->coloured().size() + clnew->antiColoured().size()) > 1 ){
        if( (clold->coloured().size() + clold->antiColoured().size()) > 1 ){
          //join the colour lines
          //I don't use the join method, because potentially only (anti)coloured
          //particles belong to one colour line
	  if(clold!=clnew){//procs are not already connected
	    while ( !clnew->coloured().empty() ) {
	      tPPtr p = clnew->coloured()[0];
	      clnew->removeColoured(p);
	      clold->addColoured(p);
	    }
	    while ( !clnew->antiColoured().empty() ) {
	      tPPtr p = clnew->antiColoured()[0];
	      clnew->removeAntiColoured(p);
	      clold->addAntiColoured(p);
	    }
	  }

        }else{
          //if pold is the only member on it's 
          //colour line, remove it.
          clold->removeColoured(pold, anti);
          //and add it to clnew
          clnew->addColoured(pold, anti);
        }    
      }else{//pnnew is the only member on it's colour line.
        clnew->removeColoured(pnew, !anti);
        clold->addColoured(pnew, !anti);
      }
    }else{//there is no coline at all for pnew
      clold->addColoured(pnew, !anti);
    }
    //end of loop
  }
  return;
}

void HwRemDecayer::doSplit(pair<tPPtr, tPPtr> partons, bool first) {
  if(!theSplittingOnOff) return;
  
  try{
    split(partons.first, theContent.first, theRems.first, 
	  theUsed.first, theMaps.first, first);
  }catch(ShowerHandler::ExtraScatterVeto){
    theX.first -= partons.first->momentum().rho()/parent(theRems.first)->momentum().rho();
    throw ShowerHandler::ExtraScatterVeto();
  }
  
  try{
    split(partons.second, theContent.second, theRems.second, 
	  theUsed.second, theMaps.second, first);
  }catch(ShowerHandler::ExtraScatterVeto){
    //case of the first forcedSplitting worked fine
    theX.first -= partons.first->momentum().rho()/parent(theRems.first)->momentum().rho();
    theX.second -= partons.second->momentum().rho()/parent(theRems.second)->momentum().rho();
    
    //case of the first interaction
    //throw veto immediately, because event get rejected anyway.
    if(first) throw ShowerHandler::ExtraScatterVeto();

    //secondary interactions have to end on a gluon, if parton 
    //was NOT a gluon, the forced splitting particles must be removed
    if(partons.first->id() != ParticleID::g){
      if(partons.first==theMaps.first.back().first) 
	theUsed.first -= theMaps.first.back().second->momentum();
      else
	theUsed.first -= theMaps.first.back().first->momentum();
      
      thestep->removeParticle(theMaps.first.back().first);
      thestep->removeParticle(theMaps.first.back().second);
    }
    theMaps.first.pop_back();
    throw ShowerHandler::ExtraScatterVeto();
  }
}

void HwRemDecayer::finalize(){
  if(!theSplittingOnOff) return;

  PPtr diquark;
  //Do the final Rem->Diquark or Rem->quark "decay"
  diquark = theForcedSplitter->
    finalSplit(theRems.first, theContent.first.RemID(), 
	       theUsed.first, thestep);
  theMaps.first.push_back(make_pair(diquark, tPPtr()));

  diquark = theForcedSplitter->
    finalSplit(theRems.second, theContent.second.RemID(), 
	       theUsed.second, thestep);
  theMaps.second.push_back(make_pair(diquark, tPPtr()));

  setRemMasses();
  fixColours(theMaps.first, theanti.first);
  fixColours(theMaps.second, theanti.second);
}

ParticleVector HwRemDecayer::decay(const DecayMode &, 
				   const Particle &, Step &) const {
  throw Exception() << "HwRemDecayer::decay(...) "
		    << "must not be called explicitely."
		    << Exception::runerror; 
  return PVector();
}


void HwRemDecayer::persistentOutput(PersistentOStream & os) const {
  os << theForcedSplitter << theSplittingOnOff;
}

void HwRemDecayer::persistentInput(PersistentIStream & is, int) {
  is >> theForcedSplitter >> theSplittingOnOff;
}

ClassDescription<HwRemDecayer> HwRemDecayer::initHwRemDecayer;
// Definition of the static class description member.

void HwRemDecayer::Init() {

  static ClassDocumentation<HwRemDecayer> documentation
    ("There is no documentation for the HwRemDecayer class");

  static Reference<HwRemDecayer,ForcedSplitting> interfaceForcedSplitting
    ("ForcedSplitting",
     "Object used for the forced splitting of the Remnant",
     &HwRemDecayer::theForcedSplitter, false, false, true, false, false);

  static Switch<HwRemDecayer,bool> interfaceSplittingOnOff
    ("SplittingOnOff", "flag to switch the ForcedSplitting on or off, i.e. switch the entire class on or off",
     &HwRemDecayer::theSplittingOnOff, 1, false, false);

  static SwitchOption interfaceSplittingOnOff0
    (interfaceSplittingOnOff,"ForcedSplitting-OFF","Forced Splitting is OFF", 0);
  static SwitchOption interfaceSplittingOnOff1
    (interfaceSplittingOnOff,"ForcedSplitting-ON","Forced Splitting is ON", 1);

}

