// -*- C++ -*-
//
// HwRemDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HwRemDecayer class.
//

#include "HwRemDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <ThePEG/Interface/Reference.h>  
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Interface/Switch.h>
#include "Herwig++/Shower/ShowerHandler.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/Throw.h"

using namespace Herwig;

void HwRemDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_kinCutoff, GeV) << _range 
     << _zbin << _ybin << _nbinmax << _alpha << DISRemnantOpt_;
}

void HwRemDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_kinCutoff, GeV) >> _range 
     >> _zbin >> _ybin >> _nbinmax >> _alpha >> DISRemnantOpt_;
}

ClassDescription<HwRemDecayer> HwRemDecayer::initHwRemDecayer;
// Definition of the static class description member.

void HwRemDecayer::Init() {

  static ClassDocumentation<HwRemDecayer> documentation
    ("The HwRemDecayer class decays the remnant for Herwig++");

  static Parameter<HwRemDecayer,double> interfaceZBinSize
    ("ZBinSize",
     "The size of the vbins in z for the interpolation of the splitting function.",
     &HwRemDecayer::_zbin, 0.05, 0.001, 0.1,
     false, false, Interface::limited);

  static Parameter<HwRemDecayer,int> interfaceMaxBin
    ("MaxBin",
     "Maximum number of z bins",
     &HwRemDecayer::_nbinmax, 100, 10, 1000,
     false, false, Interface::limited);

  static Reference<HwRemDecayer,ShowerAlpha> interfaceAlphaS
    ("AlphaS",
     "Pointer to object to calculate the strong coupling",
     &HwRemDecayer::_alpha, false, false, true, false, false);

  static Parameter<HwRemDecayer,Energy> interfaceKinCutoff
    ("KinCutoff",
     "Parameter kinCutoff used to constrain qtilde",
     &HwRemDecayer::_kinCutoff, GeV, 0.75*GeV, 0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<HwRemDecayer,double> interfaceEmissionRange
    ("EmissionRange",
     "Factor above the minimum possible value in which the forced splitting is allowed.",
     &HwRemDecayer::_range, 1.1, 1.0, 10.0,
     false, false, Interface::limited);


  static Switch<HwRemDecayer,unsigned int> interfaceDISRemnantOption
    ("DISRemnantOption",
     "Options for the treatment of the remnant in DIS",
     &HwRemDecayer::DISRemnantOpt_, 0, false, false);
  static SwitchOption interfaceDISRemnantOptionDefault
    (interfaceDISRemnantOption,
     "Default",
     "Use the minimum number of particles needed to take the recoil"
     " and allow the lepton to be used if needed",
     0);
  static SwitchOption interfaceDISRemnantOptionNoLepton
    (interfaceDISRemnantOption,
     "NoLepton",
     "Use the minimum number of particles needed to take the recoil but"
     " veto events where the lepton kinematics would need to be altered",
     1);
  static SwitchOption interfaceDISRemnantOptionAllParticles
    (interfaceDISRemnantOption,
     "AllParticles",
     "Use all particles in the colour connected system to take the recoil"
     " and use the lepton if needed.",
     2);
  static SwitchOption interfaceDISRemnantOptionAllParticlesNoLepton
    (interfaceDISRemnantOption,
     "AllParticlesNoLepton",
     "Use all the particles in the colour connected system to take the"
     " recoil but don't use the lepton.",
     3);

}

ParticleVector HwRemDecayer::decay(const DecayMode &, 
				   const Particle &, Step &) const {
  throw Exception() << "HwRemDecayer::decay(...) "
		    << "must not be called explicitely."
		    << Exception::runerror;
}

void HwRemDecayer::initialize(pair<tRemPPtr, tRemPPtr> rems, tPPair beam, Step & step,
			      Energy forcedSplitScale) {
  // the step
  thestep = &step;
  // valence content of the hadrons
  theContent.first  = getHadronContent(beam.first);
  theContent.second = getHadronContent(beam.second);
  // momentum extracted from the hadrons
  theUsed.first  = Lorentz5Momentum();
  theUsed.second = Lorentz5Momentum();
  theMaps.first.clear();
  theMaps.second.clear();
  theX.first = 0.0;
  theX.second = 0.0;
  theRems = rems;
  _forcedSplitScale = forcedSplitScale;
  // check remnants attached to the right hadrons
  if( (theRems.first  && parent(theRems.first ) != beam.first ) ||
      (theRems.second && parent(theRems.second) != beam.second) )
    throw Exception() << "Remnant order wrong in "
		      << "HwRemDecayer::initialize(...)"
		      << Exception::runerror; 
  return;
}

void HwRemDecayer::doSplit(pair<tPPtr, tPPtr> partons, pair<tcPDFPtr, tcPDFPtr> pdfs,
			   bool first) {
  if(theRems.first) {
    ParticleVector children=theRems.first->children();
    for(unsigned int ix=0;ix<children.size();++ix) {
      if(children[ix]->dataPtr()==theRems.first->dataPtr()) 
	theRems.first = dynamic_ptr_cast<RemPPtr>(children[ix]);
    }
  }
  if(theRems.second) {
    ParticleVector children=theRems.second->children();
    for(unsigned int ix=0;ix<children.size();++ix) {
      if(children[ix]->dataPtr()==theRems.second->dataPtr()) 
	theRems.second = dynamic_ptr_cast<RemPPtr>(children[ix]);
    }
  }
  // forced splitting for first parton
  if(partons.first->data().coloured()) {
    try {
      split(partons.first, theContent.first, theRems.first, 
	    theUsed.first, theMaps.first, pdfs.first, first);
    }
    catch(ShowerHandler::ExtraScatterVeto) {
      theX.first -= partons.first->momentum().rho()/
	parent(theRems.first)->momentum().rho();
      throw ShowerHandler::ExtraScatterVeto();
    }
  }
  // forced splitting for second parton
  if(partons.second->data().coloured()) {   
    try{
      split(partons.second, theContent.second, theRems.second, 
	    theUsed.second, theMaps.second, pdfs.second, first);
      // additional check for the remnants
      // if can't do the rescale veto the emission
      if(!first&&partons.first->data().coloured()&&
	 partons.second->data().coloured()) {
	Lorentz5Momentum pnew[2]=
	  {theRems.first->momentum()  - theUsed.first  - partons.first->momentum(),
	   theRems.second->momentum() - theUsed.second - partons.second->momentum()};
	
	pnew[0].setMass(getParticleData(theContent.first.RemID())->constituentMass());
	pnew[0].rescaleEnergy();
	pnew[1].setMass(getParticleData(theContent.second.RemID())->constituentMass());
	pnew[1].rescaleEnergy();
	
	for(unsigned int iy=0; iy<theRems.first->children().size(); ++iy)
	  pnew[0] += theRems.first->children()[iy]->momentum();
     
	for(unsigned int iy=0; iy<theRems.second->children().size(); ++iy)
	  pnew[1] += theRems.second->children()[iy]->momentum();

	Lorentz5Momentum ptotal=
	  theRems.first ->momentum()-partons.first ->momentum()+
	  theRems.second->momentum()-partons.second->momentum();
     
	if(ptotal.m() < (pnew[0].m() + pnew[1].m()) ) {
	  if(partons.second->id() != ParticleID::g){
	    if(partons.second==theMaps.second.back().first) 
	      theUsed.second -= theMaps.second.back().second->momentum();
	    else
	      theUsed.second -= theMaps.second.back().first->momentum();
         
	    thestep->removeParticle(theMaps.second.back().first);
	    thestep->removeParticle(theMaps.second.back().second);
	  }
	  theMaps.second.pop_back();
	  throw ShowerHandler::ExtraScatterVeto();
	}
      }
    }
    catch(ShowerHandler::ExtraScatterVeto){
      if(!partons.first||!partons.second||
	 !theRems.first||!theRems.second) 
	throw ShowerHandler::ExtraScatterVeto();
      //case of the first forcedSplitting worked fine
      theX.first -= partons.first->momentum().rho()/
	parent(theRems.first)->momentum().rho();
      theX.second -= partons.second->momentum().rho()/
	parent(theRems.second)->momentum().rho();

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
}

void HwRemDecayer::finalize() {
  PPtr diquark;
  //Do the final Rem->Diquark or Rem->quark "decay"
  if(theRems.first) {
    diquark = finalSplit(theRems.first, theContent.first.RemID(), 
			 theUsed.first);
    theMaps.first.push_back(make_pair(diquark, tPPtr()));
  }
  if(theRems.second) {
    diquark = finalSplit(theRems.second, theContent.second.RemID(), 
			 theUsed.second);
    theMaps.second.push_back(make_pair(diquark, tPPtr()));
  }
  setRemMasses();
  if(theRems.first)  fixColours(theMaps.first , theanti.first );
  if(theRems.second) fixColours(theMaps.second, theanti.second);
}

HwRemDecayer::HadronContent
HwRemDecayer::getHadronContent(tcPPtr hadron) const {
  HadronContent hc;
  hc.hadron = hadron->dataPtr();
  long id(hadron->id());
  // baryon
  if(BaryonMatcher::Check(hadron->data())) {
    hc.sign = id < 0? -1: 1;
    hc.flav.push_back((id = abs(id)/10)%10);
    hc.flav.push_back((id /= 10)%10);
    hc.flav.push_back((id /= 10)%10);
    hc.extracted = -1;
  }
  else if(hadron->data().id()==ParticleID::gamma) {
    hc.sign = 1;
    for(int ix=1;ix<6;++ix) {
      hc.flav.push_back( ix);
      hc.flav.push_back(-ix);
    }
  }
  return hc;
}

long HwRemDecayer::HadronContent::RemID() const{
  if(extracted == -1)
    throw Exception() << "Try to build a Diquark id without "
		      << "having extracted something in "
		      << "HwRemDecayer::RemID(...)"
		      << Exception::runerror;
  //the hadron was a meson or photon
  if(flav.size()==2) return sign*flav[(extracted+1)%2];
 
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

void HwRemDecayer::split(tPPtr parton, HadronContent & content, 
			 tRemPPtr rem, Lorentz5Momentum & used, 
			 PartnerMap &partners, tcPDFPtr pdf, bool first) {
  theBeam = parent(rem);
  theBeamData = dynamic_ptr_cast<Ptr<BeamParticleData>::const_pointer>
    (theBeam->dataPtr());

  if(rem==theRems.first)
    theX.first  += parton->momentum().rho()/theBeam->momentum().rho();
  else
    theX.second += parton->momentum().rho()/theBeam->momentum().rho();

  double check = rem==theRems.first ? theX.first : theX.second;

  if(1.0-check < 1e-3) throw ShowerHandler::ExtraScatterVeto();

  bool anti;
  Lorentz5Momentum lastp(parton->momentum());
  int lastID(parton->id());
  ColinePtr cl;
  Energy oldQ(_forcedSplitScale);
  _pdf = pdf;
  //do nothing if already valence quark
  if(first && content.isValenceQuark(parton)) { 
    //set the extracted value, because otherwise no RemID could be generated.
    content.extract(lastID);
    // add the particle to the colour partners
    partners.push_back(make_pair(parton, tPPtr()));
    //set the sign
    anti = parton->hasAntiColour();
    if(rem==theRems.first) theanti.first  = anti;
    else                   theanti.second = anti;
    return;
  }
  //or gluon for secondaries
  else if(!first && lastID == ParticleID::g){
    partners.push_back(make_pair(parton, tPPtr()));
    return; 
  }
  // if a sea quark.antiquark forced splitting to a gluon
  // Create the new parton with its momentum and parent/child relationship set
  PPtr newSea;
  if( lastID != ParticleID::g ) {
    newSea = forceSplit(rem, -lastID, oldQ, 
			rem==theRems.first ? theX.first : theX.second,
			lastp, used,content);
    cl = new_ptr(ColourLine());
    if(newSea->id() > 0) cl->addColoured(newSea);
    else cl->addAntiColoured(newSea);
    // if a secondard scatter finished so return
    if(!first){
      partners.push_back(make_pair(parton, newSea));
      return;
    }
  }
  // otherwise evolve back to valence
  if( !content.isValenceQuark(parton) ) {
    // final valence splitting
    PPtr newValence = forceSplit(rem, 0, oldQ,
				 rem==theRems.first ? theX.first : theX.second,
				 lastp, used, content);
    // extract from the hadron to allow remnant to be determined
    content.extract(newValence->id());
    // case of a gluon going into the hard subprocess
    if( lastID == ParticleID::g ) {
      partners.push_back(make_pair(parton, tPPtr()));
      anti = newValence->hasAntiColour();
      if(rem==theRems.first) theanti.first  = anti;
      else                   theanti.second = anti;
      parton->colourLine(!anti)->addColoured(newValence, anti);
      return;
    }
    //The valence quark will always be connected to the sea quark with opposite sign
    tcPPtr particle;
    if(lastID*newValence->id() < 0){
      particle = parton;
      partners.push_back(make_pair(newSea, tPPtr()));
    } 
    else {
      particle = newSea;
      partners.push_back(make_pair(parton, tPPtr()));
    }
    anti = newValence->hasAntiColour();
    if(rem==theRems.first) theanti.first  = anti;
    else                   theanti.second = anti;
    
    if(particle->colourLine()) 
      particle->colourLine()->addAntiColoured(newValence);
    if(particle->antiColourLine()) 
      particle->antiColourLine()->addColoured(newValence);  
  }
  return;
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

PPtr HwRemDecayer::forceSplit(const tRemPPtr rem, long child, Energy &lastQ, 
			      double &lastx, Lorentz5Momentum &pf, 
			      Lorentz5Momentum &p, HadronContent & content) const {
  // beam momentum
  Lorentz5Momentum beam = theBeam->momentum();
  // the last scale is minimum of last value and upper limit
  Energy minQ=_range*_kinCutoff*sqrt(lastx)/(1-lastx);
  if(minQ>lastQ) lastQ=minQ;
  // generate the new value of qtilde
  // weighted towards the lower value: dP/dQ = 1/Q -> Q(R) =
  // Q0 (Qmax/Q0)^R
  Energy q;
  double zmin,zmax,yy;
  do {
    q = minQ*pow(lastQ/minQ,UseRandom::rnd());
    zmin = lastx;
    yy   = 1.+0.5*sqr(_kinCutoff/q);
    zmax = yy - sqrt(sqr(yy)-1.); 
  }
  while(zmax<zmin);
  // now generate z as in FORTRAN HERWIG
  // use y = ln(z/(1-z)) as integration variable
  double ymin=log(zmin/(1.-zmin));
  double ymax=log(zmax/(1.-zmax));
  double dely=ymax-ymin;
  unsigned int nz=_nbinmax;
  dely/=nz;
  yy=ymin+0.5*dely;
  vector<int> ids;
  if(child!=0) ids.push_back(ParticleID::g);
  else         {
    ids=content.flav;
    for(unsigned int ix=0;ix<ids.size();++ix) ids[ix] *= content.sign;
  }
  // probabilities of the different types of possible splitting
  map<long,pair<double,vector<double> > > partonprob;
  double ptotal(0.);
  for(unsigned int iflav=0;iflav<ids.size();++iflav) {
    // only do each parton once
    if(partonprob.find(ids[iflav])!=partonprob.end()) continue;
    // particle data object
    tcPDPtr in = getParticleData(ids[iflav]);
    double psum(0.);
    vector<double> prob;
    for(unsigned int iz=0;iz<nz;++iz) {
      double ez=exp(yy);
      double wr=1.+ez;
      double zr=wr/ez;
      double wz=1./wr;
      double zz=wz*ez;
      double az=wz*zz*_alpha->value(sqr(max(wz*q,_kinCutoff)));
      // g -> q qbar
      if(ids[iflav]==ParticleID::g) {
	// calculate splitting function   
	// SP as q is always less than forcedSplitScale, the pdf scale is fixed 
	// pdfval = _pdf->xfx(theBeamData,in,sqr(q),lastx*zr); 
	double pdfval=_pdf->xfx(theBeamData,in,sqr(_forcedSplitScale),lastx*zr);
	psum += pdfval*az*0.5*(sqr(zz)+sqr(wz));
      }
      // q -> q g
      else {
	// calculate splitting function
	// SP as q is always less than forcedSplitScale, the pdf scale is fixed 
 	// pdfval = _pdf->xfx(theBeamData,in,sqr(q),lastx*zr);
 	double pdfval=_pdf->xfx(theBeamData,in,sqr(_forcedSplitScale),lastx*zr);
	psum += pdfval*az*4./3.*(1.+sqr(wz))*zr;
      }
      prob.push_back(psum);
      yy+=dely;
    }
    partonprob[ids[iflav]] = make_pair(psum,prob);
    ptotal+=psum;
  }
  // select the flavour
  ptotal *= UseRandom::rnd();
  map<long,pair<double,vector<double> > >::const_iterator pit;
  for(pit=partonprob.begin();pit!=partonprob.end();++pit) {
    if(pit->second.first>=ptotal) break;
    else                          ptotal -= pit->second.first;
  }
  if(pit==partonprob.end()) 
    throw Exception() << "Can't select parton for forced backward evolution in "
		      << "HwRemDecayer::forceSplit" << Exception::eventerror;
  // select z 
  unsigned int iz=0;
  for(;iz<pit->second.second.size();++iz) {
    if(pit->second.second[iz]>ptotal) break;
  }
  if(iz==pit->second.second.size()) --iz;
  double ey=exp(ymin+dely*(float(iz+1)-UseRandom::rnd()));
  double z=ey/(1.+ey);
  Energy2 pt2=sqr((1.-z)*q)- z*sqr(_kinCutoff);
  // create the particle
  if(pit->first!=ParticleID::g) child=pit->first;
  PPtr parton = getParticleData(child)->produceParticle();
  Energy2 emittedm2 = sqr(parton->dataPtr()->constituentMass());
  // Now boost pcm and pf to z only frame
  Lorentz5Momentum p_ref = Lorentz5Momentum(0.0*MeV,  beam.vect());
  Lorentz5Momentum n_ref = Lorentz5Momentum(0.0*MeV, -beam.vect());
  // generate phi and compute pt of branching
  double phi = Constants::twopi*UseRandom::rnd();
  Energy pt=sqrt(pt2);
  Lorentz5Momentum qt   = LorentzMomentum(pt*cos(phi), pt*sin(phi), 0.0*MeV, 0.0*MeV);
  // compute alpha for previous particle
  Energy2 p_dot_n  = p_ref*n_ref;
  double lastalpha =    pf*n_ref/p_dot_n;
  Lorentz5Momentum qtout=qt;
  Energy2 qtout2=-qt*qt;
  double alphaout=(1.-z)/z*lastalpha;
  double betaout=0.5*(emittedm2+qtout2)/alphaout/p_dot_n;
  Lorentz5Momentum k=alphaout*p_ref+betaout*n_ref+qtout;
  k.rescaleMass();
  parton->set5Momentum(k);
  pf+=k;
  lastQ=q;
  lastx/=z;
  p += parton->momentum();
  thestep->addDecayProduct(rem,parton,false);
  return parton;
}

void HwRemDecayer::setRemMasses() const {
  // get the masses of the remnants
  Energy mrem[2];
  Lorentz5Momentum ptotal,pnew[2];
  vector<tRemPPtr> theprocessed;
  theprocessed.push_back(theRems.first);
  theprocessed.push_back(theRems.second);
  // one remnant in e.g. DIS
  if(!theprocessed[0]||!theprocessed[1]) {
    tRemPPtr rem = theprocessed[0] ? theprocessed[0] : theprocessed[1];
    Lorentz5Momentum deltap(rem->momentum());
    // find the diquark and momentum we still need in the energy
    tPPtr diquark;
    vector<PPtr> progenitors;
    for(unsigned int ix=0;ix<rem->children().size();++ix) {
      if(!DiquarkMatcher::Check(rem->children()[ix]->data())) {
	progenitors.push_back(rem->children()[ix]);
	deltap -= rem->children()[ix]->momentum();
      }
      else 
	diquark = rem->children()[ix];
    }
    // now find the total momentum of the hadronic final-state to
    // reshuffle against
    // find the hadron for this remnant
    tPPtr hadron=rem;
    do hadron=hadron->parents()[0];
    while(!hadron->parents().empty());
    // find incoming parton to hard process from this hadron
    tPPtr hardin = 
      generator()->currentEvent()->primaryCollision()->incoming().first==hadron ?
      generator()->currentEvent()->primarySubProcess()->incoming().first :
      generator()->currentEvent()->primarySubProcess()->incoming().second;
    tPPtr parent=hardin;
    vector<PPtr> tempprog;
    // find the outgoing particles emitted from the backward shower
    do {
      assert(!parent->parents().empty());
      tPPtr newparent=parent->parents()[0];
      if(newparent==hadron) break;
      for(unsigned int ix=0;ix<newparent->children().size();++ix) {
	if(newparent->children()[ix]!=parent)
	  findChildren(newparent->children()[ix],tempprog);
      }
      parent=newparent;
    }
    while(parent!=hadron);
    // add to list of potential particles to reshuffle against in right order
    for(unsigned int ix=tempprog.size();ix>0;--ix) progenitors.push_back(tempprog[ix-1]);
    // final-state particles which are colour connected
    tColinePair lines = make_pair(hardin->colourLine(),hardin->antiColourLine());
    vector<PPtr> others;
    for(ParticleVector::const_iterator 
	  cit = generator()->currentEvent()->primarySubProcess()->outgoing().begin();
	cit!= generator()->currentEvent()->primarySubProcess()->outgoing().end();++cit) {
      // colour connected
      if(lines.first&&lines.first==(**cit).colourLine()) {
	findChildren(*cit,progenitors);
	continue;
      }
      // anticolour connected
      if(lines.second&&lines.second==(**cit).antiColourLine()) {
	findChildren(*cit,progenitors);
	continue;
      }
      // not connected
      for(unsigned int ix=0;ix<(**cit).children().size();++ix)
	others.push_back((**cit).children()[ix]);
    }
    // work out how much of the system needed for rescaling
    unsigned int iloc=0;
    Lorentz5Momentum psystem,ptotal;
    do {
      psystem+=progenitors[iloc]->momentum();
      ptotal = psystem + deltap;
      ptotal.rescaleMass();
      psystem.rescaleMass();
      ++iloc;
      if(ptotal.mass() > psystem.mass() + diquark->mass() &&
	 DISRemnantOpt_<2) break;
    }
    while(iloc<progenitors.size());
    if(ptotal.mass() > psystem.mass() + diquark->mass()) --iloc;
    if(iloc==progenitors.size()) {
      // try touching the lepton in dis as a last restort
      for(unsigned int ix=0;ix<others.size();++ix) {
	progenitors.push_back(others[ix]);
	psystem+=progenitors[iloc]->momentum();
	ptotal = psystem + deltap;
	ptotal.rescaleMass();
	psystem.rescaleMass();
	++iloc;
      }
      --iloc;
      if(ptotal.mass() > psystem.mass() + diquark->mass()) {
	if(DISRemnantOpt_==0||DISRemnantOpt_==2)
	  Throw<Exception>() << "Warning had to adjust the momentum of the"
			     << " non-colour connected"
			     << " final-state, e.g. the scattered lepton in DIS"
			     << Exception::warning;
	else
	  throw Exception() << "Can't set remnant momentum without adjusting "
			    << "the momentum of the"
			    << " non-colour connected"
			    << " final-state, e.g. the scattered lepton in DIS"
			    << " vetoing event"
			    << Exception::eventerror;
      }
      else {
	throw Exception() << "Can't put the remnant on-shell in HwRemDecayer::setRemMasses()"
			  << Exception::eventerror;
      }
    }
    psystem.rescaleMass();
    LorentzRotation R = Utilities::getBoostToCM(make_pair(psystem, deltap));
    Energy pz = SimplePhaseSpace::getMagnitude(sqr(ptotal.mass()), 
					       psystem.mass(), diquark->mass());
    LorentzRotation Rs(-(R*psystem).boostVector());
    Rs.boost(0.0, 0.0, pz/sqrt(sqr(pz) + sqr(psystem.mass())));
    Rs = Rs*R;
    // put remnant on shell
    deltap.transform(R);
    deltap.setMass(diquark->mass());
    deltap.setE(sqrt(sqr(diquark->mass())+sqr(pz)));
    deltap.rescaleRho();
    R.invert();
    deltap.transform(R);
    Rs = R*Rs;
    // apply transformation to required particles to absorb recoil
    for(unsigned int ix=0;ix<=iloc;++ix) {
      progenitors[ix]->deepTransform(Rs);
    }
    diquark->set5Momentum(deltap);
  }
  // two remnants
  else {
    for(unsigned int ix=0;ix<2;++ix) {
      if(!theprocessed[ix]) continue;
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
    if(ptotal.m()< (pnew[0].m()+pnew[1].m())) 
      throw Exception() << "Not enough energy in both remnants in " 
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
}

void HwRemDecayer::findChildren(tPPtr part,vector<PPtr> & particles) const {
  if(part->children().empty()) particles.push_back(part);
  else {
    for(unsigned int ix=0;ix<part->children().size();++ix)
      findChildren(part->children()[ix],particles);
  }
}
