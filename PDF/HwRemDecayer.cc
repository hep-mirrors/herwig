// -*- C++ -*-
//
// HwRemDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HwRemDecayer class.
//

#include "HwRemDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h"  
#include "ThePEG/Interface/Parameter.h"  
#include "ThePEG/Interface/Switch.h"  
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/Throw.h"

#include "Herwig/Shower/ShowerHandler.h"

using namespace Herwig;

namespace{

const bool dbg = false;


  void reShuffle(Lorentz5Momentum &p1, Lorentz5Momentum &p2, Energy m1, Energy m2){

    Lorentz5Momentum ptotal(p1+p2);
    ptotal.rescaleMass();

    if( ptotal.m() < m1+m2 ) {
      if(dbg)
	cerr << "Not enough energy to perform reshuffling \n";
      throw HwRemDecayer::ExtraSoftScatterVeto();
    }

    Boost boostv = -ptotal.boostVector();
    ptotal.boost(boostv);

    p1.boost(boostv);
    // set the masses and energies,
    p1.setMass(m1);
    p1.setE(0.5/ptotal.m()*(ptotal.m2()+sqr(m1)-sqr(m2)));
    p1.rescaleRho();
    // boost back to the lab
    p1.boost(-boostv);

    p2.boost(boostv);
    // set the masses and energies,
    p2.setMass(m2);
    p2.setE(0.5/ptotal.m()*(ptotal.m2()+sqr(m2)-sqr(m1)));
    p2.rescaleRho();
    // boost back to the lab
    p2.boost(-boostv);
  }
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

void HwRemDecayer::split(tPPtr parton, HadronContent & content, 
			 tRemPPtr rem, Lorentz5Momentum & used, 
			 PartnerMap &partners, tcPDFPtr pdf, bool first) {
  theBeam = parent(rem);
  theBeamData = dynamic_ptr_cast<Ptr<BeamParticleData>::const_pointer>
    (theBeam->dataPtr());
  double currentx = parton->momentum().rho()/theBeam->momentum().rho();
  double check = rem==theRems.first ? theX.first : theX.second;
  check += currentx;
  if(1.0-check < 1e-3) throw ShowerHandler::ExtraScatterVeto();
  bool anti;
  Lorentz5Momentum lastp(parton->momentum());
  int lastID(parton->id());
  Energy oldQ(_forcedSplitScale);
  _pdf = pdf;
  //do nothing if already valence quark
  if(first && content.isValenceQuark(parton)) { 
    //set the extracted value, because otherwise no RemID could be generated.
    content.extract(lastID);
    // add the particle to the colour partners
    partners.push_back(make_pair(parton, tPPtr()));
    //set the sign
    anti = parton->hasAntiColour() && parton->id()!=ParticleID::g;
    if(rem==theRems.first) theanti.first  = anti;
    else                   theanti.second = anti;
    // add the x and return
    if(rem==theRems.first) theX.first  += currentx;
    else                   theX.second += currentx;
    return;
  }
  //or gluon for secondaries
  else if(!first && lastID == ParticleID::g) {
    partners.push_back(make_pair(parton, tPPtr()));
    // add the x and return
    if(rem==theRems.first) theX.first  += currentx;
    else                   theX.second += currentx;
    return; 
  }
  // if a sea quark.antiquark forced splitting to a gluon
  // Create the new parton with its momentum and parent/child relationship set
  PPtr newSea;
  if( !(lastID == ParticleID::g || 
	lastID == ParticleID::gamma) ) {
    newSea = forceSplit(rem, -lastID, oldQ, currentx, lastp, used,content);
    ColinePtr cl = new_ptr(ColourLine());
    if(newSea->id() > 0) cl->    addColoured(newSea);
    else                 cl->addAntiColoured(newSea);
    // if a secondard scatter finished so return
    if(!first || content.isValenceQuark(ParticleID::g) ){
      partners.push_back(make_pair(parton, newSea));
      // add the x and return
      if(rem==theRems.first) theX.first  += currentx;
      else                   theX.second += currentx;
      if(first) content.extract(ParticleID::g);
      return;
    }
  }
  // otherwise evolve back to valence
  // final valence splitting
  PPtr newValence = forceSplit(rem, 
			       lastID!=ParticleID::gamma ? 
			       ParticleID::g : ParticleID::gamma,
			       oldQ, currentx , lastp, used, content);
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
  else if( lastID == ParticleID::gamma) {
    partners.push_back(make_pair(parton, newValence));
    anti = newValence->hasAntiColour();
    ColinePtr newLine(new_ptr(ColourLine()));
    newLine->addColoured(newValence, anti);
    if(rem==theRems.first) theanti.first  = anti;
    else                   theanti.second = anti;
    // add the x and return
    if(rem==theRems.first) theX.first  += currentx;
    else                   theX.second += currentx;
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
  // add the x and return
  if(rem==theRems.first) theX.first  += currentx;
  else                   theX.second += currentx;
  return;
}

void HwRemDecayer::doSplit(pair<tPPtr, tPPtr> partons,
			   pair<tcPDFPtr, tcPDFPtr> pdfs,
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
  if(isPartonic(partons.first )) { 
    try {
      split(partons.first, theContent.first, theRems.first, 
	    theUsed.first, theMaps.first, pdfs.first, first);
    }
    catch(ShowerHandler::ExtraScatterVeto) {
      throw ShowerHandler::ExtraScatterVeto();
    }
  }
  // forced splitting for second parton
  if(isPartonic(partons.second)) { 
    try {
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
	// add x limits
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
	  theX.second -= partons.second->momentum().rho()/
	    parent(theRems.second)->momentum().rho();
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
      //case of the first interaction
      //throw veto immediately, because event get rejected anyway.
      if(first) throw ShowerHandler::ExtraScatterVeto();
      //secondary interactions have to end on a gluon, if parton 
      //was NOT a gluon, the forced splitting particles must be removed
      if(partons.first->id() != ParticleID::g) {
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
  // veto if not enough energy for extraction
  if( !first &&(theRems.first ->momentum().e() - 
		partons.first ->momentum().e() < 1.0e-3*MeV ||
		theRems.second->momentum().e() - 
		partons.second->momentum().e() < 1.0e-3*MeV )) {
    if(partons.first->id() != ParticleID::g) {
      if(partons.first==theMaps.first.back().first) 
	theUsed.first -= theMaps.first.back().second->momentum();
      else
	theUsed.first -= theMaps.first.back().first->momentum();
      thestep->removeParticle(theMaps.first.back().first);
      thestep->removeParticle(theMaps.first.back().second);
    }
    theMaps.first.pop_back();
    if(partons.second->id() != ParticleID::g) {
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

void HwRemDecayer::mergeColour(tPPtr pold, tPPtr pnew, bool anti) const {
  ColinePtr clnew, clold;
  
  //save the corresponding colour lines
  clold = pold->colourLine(anti);
  clnew = pnew->colourLine(!anti);

  assert(clold);

  // There is already a colour line (not the final diquark)
  if(clnew){

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
    } else{//pnnew is the only member on it's colour line.
      clnew->removeColoured(pnew, !anti);
      clold->addColoured(pnew, !anti);
    }
  } else {//there is no coline at all for pnew
    clold->addColoured(pnew, !anti);
  }
}

void HwRemDecayer::fixColours(PartnerMap partners, bool anti, 
			      double colourDisrupt) const {
  PartnerMap::iterator prev;
  tPPtr pnew, pold;
  assert(partners.size()>=2);

  PartnerMap::iterator it=partners.begin();
  while(it != partners.end()) {
    //skip the first one to have a partner
    if(it==partners.begin()){
      it++;
      continue;
    }
    
    prev = it - 1;
    //determine the particles to work with
    pold = prev->first;
    if(prev->second) {
      if(!pold->coloured())
	pold = prev->second;
      else if(pold->hasAntiColour() != anti)
	pold = prev->second;
    }
    assert(pold);

    pnew = it->first;
    if(it->second) {
      if(it->second->colourLine(!anti)) //look for the opposite colour
	pnew = it->second;
    }
    assert(pnew);

    // Implement the disruption of colour connections
    if( it != partners.end()-1 ) {//last one is diquark-has to be connected
      
      //has to be inside the if statement, so that the probability is
      //correctly counted:
      if( UseRandom::rnd() < colourDisrupt ){

	if(!it->second){//check, whether we have a gluon
	  mergeColour(pnew, pnew, anti);
	}else{
	  if(pnew==it->first)//be careful about the order
	    mergeColour(it->second, it->first, anti);
	  else
	    mergeColour(it->first, it->second, anti);
	}

	it = partners.erase(it);
	continue;
      }

    }
    // regular merging
    mergeColour(pold, pnew, anti);
    //end of loop
    it++;
  }
  return;
}

PPtr HwRemDecayer::forceSplit(const tRemPPtr rem, long child, Energy &lastQ, 
			      double &lastx, Lorentz5Momentum &pf, 
			      Lorentz5Momentum &p, 
			      HadronContent & content) const {
  static const double eps=1e-6;
  // beam momentum
  Lorentz5Momentum beam = theBeam->momentum();
  // the last scale is minimum of last value and upper limit
  Energy minQ=_range*_kinCutoff*sqrt(lastx)/(1-lastx);
  if(minQ>lastQ) lastQ=minQ;
  // generate the new value of qtilde
  // weighted towards the lower value: dP/dQ = 1/Q -> Q(R) =
  // Q0 (Qmax/Q0)^R
  Energy q;
  unsigned int ntry=0,maxtry=100;
  double xExtracted = rem==theRems.first ? theX.first : theX.second;
  double zmin=  lastx/(1.-xExtracted) ,zmax,yy;

  if(1-lastx<eps) throw ShowerHandler::ExtraScatterVeto();
  do {
    q = minQ*pow(lastQ/minQ,UseRandom::rnd());
    yy   = 1.+0.5*sqr(_kinCutoff/q);
    zmax = yy - sqrt(sqr(yy)-1.); 
    ++ntry;
  }
  while(zmax<zmin&&ntry<maxtry);
  if(ntry==maxtry) throw ShowerHandler::ExtraScatterVeto();
  if(zmax-zmin<eps) throw ShowerHandler::ExtraScatterVeto();
  // now generate z as in FORTRAN HERWIG
  // use y = ln(z/(1-z)) as integration variable
  double ymin=log(zmin/(1.-zmin));
  double ymax=log(zmax/(1.-zmax));
  double dely=ymax-ymin;
  unsigned int nz=_nbinmax;
  dely/=nz;
  yy=ymin+0.5*dely;
  vector<int> ids;
  if(child==21||child==22) {
    ids=content.flav;
    for(unsigned int ix=0;ix<ids.size();++ix) ids[ix] *= content.sign;
  }
  else {
    ids.push_back(ParticleID::g);
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
      double coup = child!=22 ? 
	_alphaS ->value(sqr(max(wz*q,_kinCutoff))) :
	_alphaEM->value(sqr(max(wz*q,_kinCutoff)));
      double az=wz*zz*coup;
      // g -> q qbar
      if(ids[iflav]==ParticleID::g) {
	// calculate splitting function   
	// SP as q is always less than forcedSplitScale, the pdf scale is fixed 
	// pdfval = _pdf->xfx(theBeamData,in,sqr(q),lastx*zr); 
	double pdfval=_pdf->xfx(theBeamData,in,sqr(_forcedSplitScale),lastx*zr);
	if(pdfval>0.) psum += pdfval*az*0.5*(sqr(zz)+sqr(wz));
      }
      // q -> q g
      else {
	// calculate splitting function
	// SP as q is always less than forcedSplitScale, the pdf scale is fixed 
 	// pdfval = _pdf->xfx(theBeamData,in,sqr(q),lastx*zr);
 	double pdfval=_pdf->xfx(theBeamData,in,sqr(_forcedSplitScale),lastx*zr);
	if(pdfval>0.) psum += pdfval*az*4./3.*(1.+sqr(wz))*zr;
      }
      if(psum>0.) prob.push_back(psum);
      yy+=dely;
    }
    if(psum>0.) partonprob[ids[iflav]] = make_pair(psum,prob);
    ptotal+=psum;
  }
  // select the flavour
  if(ptotal==0.) throw ShowerHandler::ExtraScatterVeto();
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
  Lorentz5Momentum p_ref = Lorentz5Momentum(ZERO,  beam.vect());
  Lorentz5Momentum n_ref = Lorentz5Momentum(ZERO, -beam.vect());
  // generate phi and compute pt of branching
  double phi = Constants::twopi*UseRandom::rnd();
  Energy pt=sqrt(pt2);
  Lorentz5Momentum qt   = LorentzMomentum(pt*cos(phi), pt*sin(phi), ZERO, ZERO);
  Axis axis(p_ref.vect().unit());
  if(axis.perp2()>0.) {
    LorentzRotation rot;
    double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
    rot.setRotate(acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    qt.transform(rot);
  }
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
	 psystem.mass()>1*MeV && DISRemnantOpt_<2 && ptotal.e() > 0.*GeV ) break;
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

void HwRemDecayer::initSoftInteractions(Energy ptmin, InvEnergy2 beta){
  ptmin_ = ptmin;
  beta_ = beta;
}

Energy HwRemDecayer::softPt() const {
  Energy2 pt2(ZERO);
  double xmin(0.0), xmax(1.0), x(0);

  if(beta_ == ZERO){
    return UseRandom::rnd(0.0,(double)(ptmin_/GeV))*GeV;
  }

  if(beta_ < ZERO){
    xmin = 1.0;
    xmax = exp( -beta_*sqr(ptmin_) );    
  }else{
    xmin = exp( -beta_*sqr(ptmin_) );    
    xmax = 1.0;
  }
  x = UseRandom::rnd(xmin, xmax);
  pt2 = 1.0/beta_ * log(1/x);
  
  if( pt2 < ZERO || pt2 > sqr(ptmin_) )
    throw Exception() << "HwRemDecayer::softPt generation of pt "
                      << "outside allowed range [0," << ptmin_/GeV << "]."
                      << Exception::runerror; 

  return sqrt(pt2);
}

void HwRemDecayer::softKinematics(Lorentz5Momentum &r1, Lorentz5Momentum &r2, 
				  Lorentz5Momentum &g1, Lorentz5Momentum &g2) const {

  g1 = Lorentz5Momentum();
  g2 = Lorentz5Momentum();
  //All necessary variables for the two soft gluons
  Energy pt(softPt()), pz(ZERO);
  Energy2 pz2(ZERO);
  double phi(UseRandom::rnd(2.*Constants::pi));
  double x_g1(0.0), x_g2(0.0);
  
  //Get the external momenta
  tcPPair beam(generator()->currentEventHandler()->currentCollision()->incoming());
  Lorentz5Momentum P1(beam.first->momentum()), P2(beam.second->momentum());

  if(dbg){
    cerr << "new event --------------------\n"
	 << *(beam.first) << *(softRems_.first) 
	 << "-------------------\n"
	 << *(beam.second) << *(softRems_.second) << endl;
  }

  //Get x_g1 and x_g2
  //first limits
  double xmin  = sqr(ptmin_)/4.0/(P1+P2).m2();
  double x1max = (r1.e()+abs(r1.z()))/(P1.e() + abs(P1.z()));
  double x2max = (r2.e()+abs(r2.z()))/(P2.e() + abs(P2.z()));
  //now generate according to 1/x
  x_g1 = xmin * exp(UseRandom::rnd(log(x1max/xmin)));
  x_g2 = xmin * exp(UseRandom::rnd(log(x2max/xmin)));

  if(dbg)
    cerr << x1max << " " << x_g1 << endl << x2max << " " << x_g2 << endl;

  Lorentz5Momentum ig1, ig2, cmf;

  ig1 = x_g1*P1;
  ig2 = x_g2*P2;

  ig1.setMass(mg_);
  ig2.setMass(mg_);
  ig1.rescaleEnergy();
  ig2.rescaleEnergy();

  cmf = ig1 + ig2;
  //boost vector from cmf to lab
  Boost boostv(cmf.boostVector());

  //outgoing gluons in cmf
  g1.setMass(mg_);
  g2.setMass(mg_);

  g1.setX(pt*cos(phi));
  g2.setX(-pt*cos(phi));
  g1.setY(pt*sin(phi));
  g2.setY(-pt*sin(phi));
  
  pz2 = cmf.m2()/4 - sqr(mg_) - sqr(pt);

  if(pz2/GeV2 < 0.0){
    if(dbg)
      cerr << "EXCEPTION not enough energy...." << endl;
    throw ExtraSoftScatterVeto();
  }

  if(UseRandom::rndbool())
    pz = sqrt(pz2);
  else
    pz = -sqrt(pz2);

  if(dbg)
    cerr << "pz has been calculated to: " << pz/GeV << endl;


  g1.setZ(pz);
  g2.setZ(-pz);
  g1.rescaleEnergy();
  g2.rescaleEnergy();

  if(dbg){
    cerr << "check inv mass in cmf frame: " << (g1+g2).m()/GeV 
	 << " vs. lab frame: " << (ig1+ig2).m()/GeV << endl;
  }

  g1.boost(boostv);
  g2.boost(boostv);

  //recalc the remnant momenta
  Lorentz5Momentum r1old(r1), r2old(r2);
  r1 -= ig1;
  r2 -= ig2;

  try{
    reShuffle(r1, r2, r1old.m(), r2old.m());
  }catch(ExtraSoftScatterVeto){
    r1 = r1old;
    r2 = r2old;
    throw ExtraSoftScatterVeto();
  }


  if(dbg){
    cerr << "remnant 1,2 momenta: " << r1/GeV << "--" << r2/GeV << endl;
    cerr << "remnant 1,2 masses: " << r1.m()/GeV << " " << r2.m()/GeV << endl;
    cerr << "check momenta in the lab..." << (-r1old-r2old+r1+r2+g1+g2)/GeV << endl;
  }
}

void HwRemDecayer::doSoftInteractions(unsigned int N) {
  if(N == 0) return;
  if(!softRems_.first || !softRems_.second)
    throw Exception() << "HwRemDecayer::doSoftInteractions: no "
                      << "Remnants available."
                      << Exception::runerror; 

  if( ptmin_ == -1.*GeV )
    throw Exception() << "HwRemDecayer::doSoftInteractions: init "
                      << "code has not been called! call initSoftInteractions."
                      << Exception::runerror; 

  Lorentz5Momentum g1, g2;
  Lorentz5Momentum r1(softRems_.first->momentum()), r2(softRems_.second->momentum());
  unsigned int tries(1), i(0);

  for(i=0; i<N; i++){
    //check how often this scattering has been regenerated
    if(tries > maxtrySoft_) break;

    if(dbg){
      cerr << "new try \n" << *softRems_.first << *softRems_.second << endl;
    }

    try{
      softKinematics(r1, r2, g1, g2);

    }catch(ExtraSoftScatterVeto){
      tries++;
      i--;
      continue;
    }
  
    PPair oldrems = softRems_;
    PPair gluons = make_pair(addParticle(softRems_.first, ParticleID::g, g1), 
			     addParticle(softRems_.second, ParticleID::g, g2));
    //now reset the remnants with the new ones
    softRems_.first = addParticle(softRems_.first, softRems_.first->id(), r1);
    softRems_.second = addParticle(softRems_.second, softRems_.second->id(), r2);

    //do the colour connections
    pair<bool, bool> anti = make_pair(oldrems.first->hasAntiColour(),
				      oldrems.second->hasAntiColour());

    ColinePtr cl1 = new_ptr(ColourLine());
    ColinePtr cl2 = new_ptr(ColourLine());
    if( UseRandom::rnd() < colourDisrupt_ ){//this is the member variable, i.e. SOFT colour disruption
      //connect the remnants independent of the gluons
      oldrems.first->colourLine(anti.first)->addColoured(softRems_.first, anti.first);
      oldrems.second->colourLine(anti.second)->addColoured(softRems_.second, anti.second);
      //connect the gluons to each other
      cl1->addColoured(gluons.first);
      cl1->addAntiColoured(gluons.second);
      cl2->addColoured(gluons.second);
      cl2->addAntiColoured(gluons.first);	
    }else{
      //connect the remnants to the gluons
      oldrems.first->colourLine(anti.first)->addColoured(gluons.first, anti.first);
      oldrems.second->colourLine(anti.second)->addColoured(gluons.second, anti.second);
      //and the remaining colour line to the final remnant
      cl1->addColoured(softRems_.first, anti.first);
      cl1->addColoured(gluons.first, !anti.first);
      cl2->addColoured(softRems_.second, anti.second);
      cl2->addColoured(gluons.second, !anti.second);
    }

    //reset counter
    tries = 1;
  }

  if(dbg)
    cerr << "generated " << i << "th soft scatters\n";

}

void HwRemDecayer::finalize(double colourDisrupt, unsigned int softInt){
  PPair diquarks;
  //Do the final Rem->Diquark or Rem->quark "decay"
  if(theRems.first) {
    diquarks.first = finalSplit(theRems.first, theContent.first.RemID(), 
				theUsed.first);
    theMaps.first.push_back(make_pair(diquarks.first, tPPtr()));
  }
  if(theRems.second) {
    diquarks.second = finalSplit(theRems.second, theContent.second.RemID(), 
				 theUsed.second);
    theMaps.second.push_back(make_pair(diquarks.second, tPPtr()));
  }
  setRemMasses();
  if(theRems.first) {
    fixColours(theMaps.first, theanti.first, colourDisrupt);
    if(theContent.first.hadron->id()==ParticleID::pomeron&&
       pomeronStructure_==0) fixColours(theMaps.first, !theanti.first, colourDisrupt);
  }
  if(theRems.second) {
    fixColours(theMaps.second, theanti.second, colourDisrupt);
    if(theContent.second.hadron->id()==ParticleID::pomeron&&
       pomeronStructure_==0) fixColours(theMaps.second, !theanti.second, colourDisrupt);
  }

  if( !theRems.first || !theRems.second ) return;
  //stop here if we don't have two remnants
  softRems_ = diquarks;
  doSoftInteractions(softInt);
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
  else if(hadron->data().id()==ParticleID::gamma ||
	  (hadron->data().id()==ParticleID::pomeron && pomeronStructure_==1)) {
    hc.sign = 1;
    for(int ix=1;ix<6;++ix) {
      hc.flav.push_back( ix);
      hc.flav.push_back(-ix);
    }
  }
  else if(hadron->data().id()==ParticleID::pomeron ) {
    hc.sign = 1;
    hc.flav.push_back(ParticleID::g);
    hc.flav.push_back(ParticleID::g);
  }
  else if(hadron->data().id()==ParticleID::reggeon ) {
    hc.sign = 1;
    for(int ix=1;ix<3;++ix) {
      hc.flav.push_back( ix);
      hc.flav.push_back(-ix);
    }
  }
  hc.pomeronStructure = pomeronStructure_;
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


tPPtr HwRemDecayer::addParticle(tcPPtr parent, long id, Lorentz5Momentum p) const {

  PPtr newp = new_ptr(Particle(getParticleData(id)));
  newp->set5Momentum(p);
  // Add the new remnant to the step, but don't do colour connections
  thestep->addDecayProduct(parent,newp,false);
  return newp;
}

void HwRemDecayer::findChildren(tPPtr part,vector<PPtr> & particles) const {
  if(part->children().empty()) particles.push_back(part);
  else {
    for(unsigned int ix=0;ix<part->children().size();++ix)
      findChildren(part->children()[ix],particles);
  }
}

ParticleVector HwRemDecayer::decay(const DecayMode &, 
				   const Particle &, Step &) const {
  throw Exception() << "HwRemDecayer::decay(...) "
		    << "must not be called explicitely."
		    << Exception::runerror;
}

void HwRemDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_kinCutoff, GeV) << _range << _zbin << _ybin 
     << _nbinmax << _alphaS << _alphaEM << DISRemnantOpt_
     << maxtrySoft_ << colourDisrupt_ << pomeronStructure_
     << ounit(mg_,GeV) << ounit(ptmin_,GeV) << ounit(beta_,sqr(InvGeV))
     << allowTop_;
}

void HwRemDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_kinCutoff, GeV) >> _range >> _zbin >> _ybin 
     >> _nbinmax >> _alphaS >> _alphaEM >> DISRemnantOpt_
     >> maxtrySoft_ >> colourDisrupt_ >> pomeronStructure_
     >> iunit(mg_,GeV) >> iunit(ptmin_,GeV) >> iunit(beta_,sqr(InvGeV))
     >> allowTop_;
}

ClassDescription<HwRemDecayer> HwRemDecayer::initHwRemDecayer;
// Definition of the static class description member.

void HwRemDecayer::Init() {

  static ClassDocumentation<HwRemDecayer> documentation
    ("The HwRemDecayer class decays the remnant for Herwig");

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
     &HwRemDecayer::_alphaS, false, false, true, false, false);

  static Reference<HwRemDecayer,ShowerAlpha> interfaceAlphaEM
    ("AlphaEM",
     "Pointer to object to calculate the electromagnetic coupling",
     &HwRemDecayer::_alphaEM, false, false, true, false, false);

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

  static Parameter<HwRemDecayer,unsigned int> interfaceMaxTrySoft
    ("MaxTrySoft",
     "The maximum number of regeneration attempts for an additional soft scattering",
     &HwRemDecayer::maxtrySoft_, 10, 0, 100,
     false, false, Interface::limited);

  static Parameter<HwRemDecayer,double> interfacecolourDisrupt
    ("colourDisrupt",
     "Fraction of connections to additional soft subprocesses, which are colour disrupted.",
     &HwRemDecayer::colourDisrupt_, 
     1.0, 0.0, 1.0, 
     false, false, Interface::limited);


  static Switch<HwRemDecayer,unsigned int> interfacePomeronStructure
    ("PomeronStructure",
     "Option for the treatment of the valance structure of the pomeron",
     &HwRemDecayer::pomeronStructure_, 0, false, false);
  static SwitchOption interfacePomeronStructureGluon
    (interfacePomeronStructure,
     "Gluon",
     "Assume the pomeron is a two gluon state",
     0);
  static SwitchOption interfacePomeronStructureQQBar
    (interfacePomeronStructure,
     "QQBar",
     "Assumne the pomeron is q qbar as for the photon,"
     " this option is not recommended and is provide for compatiblity with POMWIG",
     1);

  static Switch<HwRemDecayer,bool> interfaceAllowTop
    ("AllowTop",
     "Allow top quarks in the hadron",
     &HwRemDecayer::allowTop_, false, false, false);
  static SwitchOption interfaceAllowTopNo
    (interfaceAllowTop,
     "No",
     "Don't allow them",
     false);
  static SwitchOption interfaceAllowTopYes
    (interfaceAllowTop,
     "Yes",
     "Allow them",
     true);

}

bool HwRemDecayer::canHandle(tcPDPtr particle, tcPDPtr parton) const {
  if(! (StandardQCDPartonMatcher::Check(*parton) || parton->id()==ParticleID::gamma) ) {
    if(abs(parton->id())==ParticleID::t) {
      if(!allowTop_)
	throw Exception() << "Top is not allow as a parton in the remant handling, please "
			  << "use a PDF which does not contain top for the remnant"
			  << " handling (preferred) or allow top in the remnant using\n"
			  << " set " << fullName() << ":AllowTop Yes\n"
			  << Exception::runerror;
    }
    else
      return false;
  }
  return HadronMatcher::Check(*particle) || particle->id()==ParticleID::gamma 
    || particle->id()==ParticleID::pomeron || particle->id()==ParticleID::reggeon;
}

bool HwRemDecayer::isPartonic(tPPtr parton) const {
  if(parton->parents().empty()) return false;
  tPPtr parent = parton->parents()[0];
  bool partonic = false;
  for(unsigned int ix=0;ix<parent->children().size();++ix) {
    if(dynamic_ptr_cast<tRemPPtr>(parent->children()[ix])) {
      partonic = true;
      break;
    }
  }
  return partonic;
}
