// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GGtoHHardGenerator class.
//

#include "GGtoHHardGenerator.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/HardTree.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Utilities/Maths.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "ThePEG/Repository/EventGenerator.h"

using namespace Herwig;

const complex<Energy2> 
GGtoHHardGenerator::_epsi = complex<Energy2>(ZERO,-1.e-10*GeV2);

GGtoHHardGenerator::GGtoHHardGenerator() : _power(2.0),_pregg(5.),
					   _preqg(3.),_pregqbar(3.),
					   _min_pt(2.*GeV),
					   _minloop(6),_maxloop(6),_massopt(0),
					   scaleFact_(1.), alphaScale_(0)
{}

void GGtoHHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << _alphaS << _power << _pregg << _preqg << _pregqbar 
     << ounit( _min_pt,GeV ) << _minloop << _maxloop << _massopt << scaleFact_
     << alphaScale_;
}

void GGtoHHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _alphaS >> _power >> _pregg >> _preqg >> _pregqbar 
     >> iunit( _min_pt, GeV )>> _minloop >> _maxloop >> _massopt >> scaleFact_
     >> alphaScale_;
}

ClassDescription<GGtoHHardGenerator> GGtoHHardGenerator::initGGtoHHardGenerator;
// Definition of the static class description member.

void GGtoHHardGenerator::Init() {

  static ClassDocumentation<GGtoHHardGenerator> documentation
    ("The GGtoHHardGenerator class implements the generation of hard QCD radiation in "
     "gg to h0 processes in the POWHEG scheme",
     "Hard QCD radiation for $gg\\to h^0$ processes in the POWHEG scheme \\cite{Hamilton:2009za}.",
     "%\\cite{Hamilton:2009za}\n"
     "\\bibitem{Hamilton:2009za}\n"
     "  K.~Hamilton, P.~Richardson and J.~Tully,\n"
     "  ``A Positive-Weight Next-to-Leading Order Monte Carlo Simulation for Higgs\n"
     "  Boson Production,''\n"
     "  JHEP {\\bf 0904}, 116 (2009)\n"
     "  [arXiv:0903.4345 [hep-ph]].\n"
     "  %%CITATION = JHEPA,0904,116;%%\n"
     );

  static Reference<GGtoHHardGenerator,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &GGtoHHardGenerator::_alphaS, false, false, true, false, false);

  static Parameter<GGtoHHardGenerator,double> interfacePower
    ("Power",
     "The power for the sampling of the matrix elements",
     &GGtoHHardGenerator::_power, 2.0, 1.0, 10.0,
     false, false, Interface::limited);

  static Parameter<GGtoHHardGenerator,double> interfacePrefactorgg
    ("Prefactorgg",
     "The prefactor for the sampling of the q qbar channel",
     &GGtoHHardGenerator::_pregg, 5.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<GGtoHHardGenerator,double> interfacePrefactorqg
    ("Prefactorqg",
     "The prefactor for the sampling of the q g channel",
     &GGtoHHardGenerator::_preqg, 3.0, 0.0, 1000.0,
     false, false, Interface::limited);
  
  static Parameter<GGtoHHardGenerator,double> interfacePrefactorgqbar
    ("Prefactorgqbar",
     "The prefactor for the sampling of the g qbar channel",
     &GGtoHHardGenerator::_pregqbar, 3.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<GGtoHHardGenerator, Energy> interfacePtMin
    ("minPt",
     "The pt cut on hardest emision generation"
     "2*(1-Beta)*exp(-sqr(intrinsicpT/RMS))/sqr(RMS)",
     &GGtoHHardGenerator::_min_pt, GeV, 2.*GeV, ZERO, 100000.0*GeV,
     false, false, Interface::limited);
  static Parameter<GGtoHHardGenerator,int> interfaceMaximumInLoop
    ("MaximumInLoop",
     "The maximum flavour of the quarks to include in the loops",
     &GGtoHHardGenerator::_maxloop, 6, 5, 6,
     false, false, Interface::limited);

  static Switch<GGtoHHardGenerator,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the treatment of the masses in the loop diagrams",
     &GGtoHHardGenerator::_massopt, 0, false, false);
  static SwitchOption interfaceMassOptionFull
    (interfaceMassOption,
     "Full",
     "Include the full mass dependence",
     0);
  static SwitchOption interfaceMassOptionLarge
    (interfaceMassOption,
     "Large",
     "Use the heavy mass limit",
     1);
  static Parameter<GGtoHHardGenerator, double> interfaceScaleFactor
    ("ScaleFactor",
     "The factor used before sHat if using a running scale",
     &GGtoHHardGenerator::scaleFact_, 1.0, 0.0, 10.0, 
     false, false, Interface::limited);
  static Switch<GGtoHHardGenerator,unsigned int> interfaceAlphaScale
    ("AlphaScale",
     "Option to use either pT or mT of the Higgs in the scale of alphaS",
     &GGtoHHardGenerator::alphaScale_, 0, false, false);
  static SwitchOption interfaceAlphaScalePT
    (interfaceAlphaScale,
     "pT",
     "Use the pt of the Higgs",
     0);
  static SwitchOption interfaceAlphaScaleMT
    (interfaceAlphaScale,
     "mT",
     "Use the mT of the Higgs",
     1);
}
   
bool GGtoHHardGenerator::canHandle(ShowerTreePtr tree) {
  // two incoming particles
  if(tree->incomingLines().size()!=2) return false;
  // should be gluons
  vector<ShowerParticlePtr> part;
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    part.push_back(cit->first->progenitor());
  }
  // check gluons
  if(part.size()!=2||part[0]->id()!=ParticleID::g||part[1]->id()!=ParticleID::g)
    return false;
  // one outgoing particles
  if(tree->outgoingLines().size()>1) return false;
  // find the outgoing particles
  part.clear();
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  for(cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    part.push_back(cjt->first->progenitor());
  }
  // check its the Higgs
  if(part[0]->id()!=ParticleID::h0) return false;
  // can handle it
  return true;
}

void GGtoHHardGenerator::doinitrun() {
  // insert the different prefactors in the vector for easy look up
  _prefactor.push_back(_pregg);
  _prefactor.push_back(_preqg);
  _prefactor.push_back(_preqg);
  _prefactor.push_back(_pregqbar);
  _prefactor.push_back(_pregqbar);
  HardestEmissionGenerator::doinitrun();
}

HardTreePtr GGtoHHardGenerator::generateHardest(ShowerTreePtr tree) {
  useMe();
  // get the particles to be showered
  _beams.clear();
  _partons.clear();
  // find the incoming particles
  ShowerParticleVector incoming;
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  vector<ShowerProgenitorPtr> particlesToShower;
  for( cit = tree->incomingLines().begin();
       cit != tree->incomingLines().end(); ++cit ) {
    incoming.push_back( cit->first->progenitor() );
    _beams.push_back( cit->first->beam() );
    _partons.push_back( cit->first->progenitor()->dataPtr() );
    particlesToShower.push_back( cit->first );
  }
  // find the higgs boson
  PPtr higgs;
  if(tree->outgoingLines().size() == 1) {
    higgs = tree->outgoingLines().begin()->first->copy();
  }
  else {
    higgs = tree->outgoingLines().begin()->first->copy()->parents()[0];
  }
  // calculate the rapidity of the higgs
  _yh = 0.5 * log((higgs->momentum().e()+higgs->momentum().z())/
 	          (higgs->momentum().e()-higgs->momentum().z()));
  _mass=higgs->mass();
  _mh2 = sqr(_mass);
  vector<Lorentz5Momentum> pnew;
  int emission_type(-1);
  // generate the hard emission and return if no emission
  if(!getEvent(pnew,emission_type)) return HardTreePtr();
  // construct the HardTree object needed to perform the showers
  ShowerParticleVector newparticles(4);
  // create the partons
  int iemit=-1;
  // g g -> h g
  if(emission_type==0) {
    newparticles[0] = new_ptr(ShowerParticle(_partons[0]      ,false));
    newparticles[1] = new_ptr(ShowerParticle(_partons[1]      ,false));
    iemit = pnew[0].z()/pnew[3].z()>0. ? 0 : 1;
  }
  // g q -> H q
  else if(emission_type==1) {
    newparticles[0] = new_ptr(ShowerParticle(_partons[0]      ,false));
    newparticles[1] = new_ptr(ShowerParticle(_out             ,false));
    iemit = 1;
  }
  // q g -> H q
  else if(emission_type==2) {
    newparticles[0] = new_ptr(ShowerParticle(_out             ,false));
    newparticles[1] = new_ptr(ShowerParticle(_partons[1]      ,false));
    iemit = 0;
  }
  // g qbar -> H qbar
  else if(emission_type==3) {
    newparticles[0] = new_ptr(ShowerParticle(_partons[0]      ,false));
    newparticles[1] = new_ptr(ShowerParticle(_out             ,false));
    iemit = 1;
  }
  // qbar g -> H qbar
  else if(emission_type==4) {
    newparticles[0] = new_ptr(ShowerParticle(_out             ,false));
    newparticles[1] = new_ptr(ShowerParticle(_partons[1]      ,false));
    iemit = 0;
  }
  // create the jet
  newparticles[3] = new_ptr(ShowerParticle(_out             , true));
  // create the boson
  newparticles[2] = new_ptr(ShowerParticle(higgs->dataPtr(),true));
  // set the momenta
  for(unsigned int ix=0;ix<4;++ix) newparticles[ix]->set5Momentum(pnew[ix]);
  // create the off-shell particle
  Lorentz5Momentum poff=pnew[iemit]-pnew[3];
  poff.rescaleMass();
  newparticles.push_back(new_ptr(ShowerParticle(_partons[iemit],false)));
  newparticles.back()->set5Momentum(poff);
  // find the sudakov for the branching
  BranchingList branchings=evolver()->splittingGenerator()->initialStateBranchings();
  long index = abs(_partons[iemit]->id());
  IdList br(3);
  // types of particle in the branching
  br[0]=newparticles[iemit]->id();
  br[1]=newparticles[  4  ]->id();
  br[2]=newparticles[  3  ]->id();
  SudakovPtr sudakov;
  for(BranchingList::const_iterator cit = branchings.lower_bound(index); 
      cit != branchings.upper_bound(index); ++cit ) {
    IdList ids = cit->second.second;
    if(ids[0]==br[0]&&ids[1]==br[1]&&ids[2]==br[2]) {
      sudakov=cit->second.first;
      break;
    }
  }
  if(!sudakov) throw Exception() << "Can't find Sudakov for the hard emission in "
 				 << "GGtoHHardGenerator::generateHardest()" 
 				 << Exception::runerror;
  vector<HardBranchingPtr> nasonin,nasonhard;
  // create the branchings for the incoming particles
  nasonin.push_back(new_ptr(HardBranching(newparticles[0],
					  iemit==0 ? sudakov : SudakovPtr(),
					  HardBranchingPtr(),HardBranching::Incoming)));
  nasonin.push_back(new_ptr(HardBranching(newparticles[1],
					  iemit==1 ? sudakov : SudakovPtr(),
					  HardBranchingPtr(),HardBranching::Incoming)));
  // create the branching for the emitted jet
  nasonin[iemit]->addChild(new_ptr(HardBranching(newparticles[3],SudakovPtr(),
						 nasonin[iemit],HardBranching::Outgoing)));
  // intermediate IS particle
  nasonhard.push_back(new_ptr(HardBranching(newparticles[4],SudakovPtr(),
					    nasonin[iemit],HardBranching::Incoming)));
  nasonin[iemit]->addChild(nasonhard.back());
  // set the colour partners
  nasonhard.back()->colourPartner(nasonin[iemit==0 ? 1 : 0]);
  nasonin[iemit==0 ? 1 : 0]->colourPartner(nasonhard.back());
  // add other particle
  nasonhard.push_back(nasonin[iemit==0 ? 1 : 0]);
  // outgoing Higgs boson
  nasonhard.push_back(new_ptr(HardBranching(newparticles[2],SudakovPtr(),
					    HardBranchingPtr(),HardBranching::Outgoing)));
  // make the tree
  HardTreePtr nasontree=new_ptr(HardTree(nasonhard,nasonin,ShowerInteraction::QCD));
  // connect the ShowerParticles with the branchings
  // and set the maximum pt for the radiation
  set<HardBranchingPtr> hard=nasontree->branchings();
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    if( _pt < _min_pt ) particlesToShower[ix]->maximumpT(_min_pt);
    else particlesToShower[ix]->maximumpT(_pt);
    for(set<HardBranchingPtr>::const_iterator mit=hard.begin();
 	mit!=hard.end();++mit) {
      if(particlesToShower[ix]->progenitor()->id()==(*mit)->branchingParticle()->id()&&
 	 (( (*mit)->status()==HardBranching::Incoming &&
	    !particlesToShower[ix]->progenitor()->isFinalState())||
	  ( (*mit)->status()==HardBranching::Outgoing&&
	    particlesToShower[ix]->progenitor()->isFinalState()))) {
	if(particlesToShower[ix]->progenitor()->momentum().z()/
	   (*mit)->branchingParticle()->momentum().z()<0.) continue;
 	nasontree->connect(particlesToShower[ix]->progenitor(),*mit);
 	if((*mit)->status()==HardBranching::Incoming) {
 	  (*mit)->beam(particlesToShower[ix]->original()->parents()[0]);
	}
 	HardBranchingPtr parent=(*mit)->parent();
 	while(parent) {
 	  parent->beam(particlesToShower[ix]->original()->parents()[0]);
 	  parent=parent->parent();
 	};
      }
    }
  }
  // calculate the shower variables
  evolver()->showerModel()->kinematicsReconstructor()->
    deconstructHardJets(nasontree,evolver(),ShowerInteraction::QCD);
  // return the answer
  return nasontree;
}

bool GGtoHHardGenerator::getEvent(vector<Lorentz5Momentum> & pnew, 
				     int & emis_type){
  // maximum pt (half of centre-of-mass energy)
  Energy maxp = 0.5*generator()->maximumCMEnergy();
  // set pt of emission to zero
  _pt=ZERO;
  //Working Variables
  Energy pt;
  double yj;
  // limits on the rapidity of the jet
  double minyj = -8.0,maxyj = 8.0;
  bool reject;
  double wgt;
  emis_type=-1;
  tcPDPtr outParton;
  for(int j=0;j<5;++j) {     
    pt = maxp;
    do {
      double a = _alphaS->overestimateValue()*_prefactor[j]*(maxyj-minyj)/(_power-1.);
      // generate next pt
      pt=GeV/pow(pow(GeV/pt,_power-1)-log(UseRandom::rnd())/a,1./(_power-1.));
      // generate rapidity of the jet
      yj=UseRandom::rnd()*(maxyj-minyj)+ minyj;
      // calculate rejection weight
      wgt=getResult(j,pt,yj,outParton);
      wgt/= _prefactor[j]*pow(GeV/pt,_power);
      reject = UseRandom::rnd()>wgt;
      //no emission event if p goes past p min - basically set to outside
      //of the histogram bounds (hopefully hist object just ignores it)
      if(pt<_min_pt){
	pt=ZERO;
	reject = false;
      }
      if(wgt>1.0) {
	ostringstream s;
	s << "GGtoHHardGenerator::getEvent weight for channel " << j
	  << "is " << wgt << " which is greater than 1";
	generator()->logWarning( Exception(s.str(), Exception::warning) );
      }
    }
    while(reject);
    // set pt of emission etc
    if(pt>_pt){
      emis_type = j;
      _pt=pt;
      _yj=yj;
      _out = outParton;
    }
  }
  //was this an (overall) no emission event?
  if(_pt<_min_pt){ 
    _pt=ZERO;
    emis_type = 5;
  }
  if(emis_type==5) return false;
  // generate the momenta of the particles
  // hadron-hadron cmf
  Energy2 s=sqr(generator()->maximumCMEnergy());
  // transverse energy
  Energy et=sqrt(_mh2+sqr(_pt));
  // first calculate all the kinematic variables
  // longitudinal real correction fractions
  double x  = _pt*exp( _yj)/sqrt(s)+et*exp( _yh)/sqrt(s);
  double y  = _pt*exp(-_yj)/sqrt(s)+et*exp(-_yh)/sqrt(s);
  // that and uhat
  // Energy2 th = -sqrt(s)*x*_pt*exp(-_yj);
  // Energy2 uh = -sqrt(s)*y*_pt*exp( _yj);
  // Energy2 sh = x*y*s;
  // reconstruct the momenta
  // incoming momenta
  pnew.push_back(Lorentz5Momentum(ZERO,ZERO,
				   x*0.5*sqrt(s), x*0.5*sqrt(s),ZERO));
  pnew.push_back(Lorentz5Momentum(ZERO,ZERO,
				  -y*0.5*sqrt(s), y*0.5*sqrt(s),ZERO));
  // outgoing momenta
  double phi(Constants::twopi*UseRandom::rnd());
  double sphi(sin(phi)),cphi(cos(phi));
  pnew.push_back(Lorentz5Momentum( cphi*_pt, sphi*_pt, et*sinh(_yh),
				   et*cosh(_yh), _mass));
  pnew.push_back(Lorentz5Momentum(-cphi*_pt,-sphi*_pt,_pt*sinh(_yj),
				  _pt*cosh(_yj),ZERO));
  return true;
}

Complex GGtoHHardGenerator::B(Energy2 s,Energy2 mf2) const {
  Complex output,pii(0.,Constants::pi);
  double rat=s/(4.*mf2);
  if(s<ZERO)
    output=2.-2.*sqrt(1.-1./rat)*log(sqrt(-rat)+sqrt(1.-rat));
  else if(s>=ZERO&&rat<1.)
    output=2.-2.*sqrt(1./rat-1.)*asin(sqrt(rat));
  else
    output=2.-sqrt(1.-1./rat)*(2.*log(sqrt(rat)+sqrt(rat-1.))-pii);
  return output;
}

complex<InvEnergy2> GGtoHHardGenerator::C(Energy2 s,Energy2 mf2) const {
  complex<InvEnergy2> output;
  Complex pii(0.,Constants::pi);
  double rat=s/(4.*mf2);
  if(s<ZERO)
    output=2.*sqr(log(sqrt(-rat)+sqrt(1.-rat)))/s;
  else if(s>=ZERO&&rat<1.)
    output=-2.*sqr(asin(sqrt(rat)))/s;
  else {
    double cosh=log(sqrt(rat)+sqrt(rat-1.));
    output=2.*(sqr(cosh)-sqr(Constants::pi)/4.-pii*cosh)/s;
  }
  return output;
}
  
Complex GGtoHHardGenerator::dIntegral(Energy2 a, Energy2 b, double y0) const {
  Complex output;
  if(b==ZERO) output=0.;
  else {
    Complex y1=0.5*(1.+sqrt(1.-4.*(a+_epsi)/b));
    Complex y2=1.-y1;
    Complex z1=y0/(y0-y1);
    Complex z2=(y0-1.)/(y0-y1);
    Complex z3=y0/(y0-y2);
    Complex z4=(y0-1.)/(y0-y2);
    output=Math::Li2(z1)-Math::Li2(z2)+Math::Li2(z3)-Math::Li2(z4);
  }
  return output;
}

complex<InvEnergy4> GGtoHHardGenerator::D(Energy2 s,Energy2 t, Energy2,
					 Energy2 mf2) const {
  Complex output,pii(0.,Constants::pi);
  Energy4 st=s*t;
  Energy4 root=sqrt(sqr(st)-4.*st*mf2*(s+t-_mh2));
  double xp=0.5*(st+root)/st,xm=1-xp;
  output = 2.*(-dIntegral(mf2,s,xp)-dIntegral(mf2,t,xp)
	       +dIntegral(mf2,_mh2,xp)+log(-xm/xp)
	       *(log((mf2+_epsi)/GeV2)-log((mf2+_epsi-s*xp*xm)/GeV2)
		 +log((mf2+_epsi-_mh2*xp*xm)/GeV2)-log((mf2+_epsi-t*xp*xm)/GeV2)));
  return output/root;
}
  
complex<Energy> GGtoHHardGenerator::me1(Energy2 s,Energy2 t,Energy2 u, Energy2 mf2,
				       unsigned int i ,unsigned int j ,unsigned int k ,
				       unsigned int i1,unsigned int j1,unsigned int k1) const {
  Energy2 s1(s-_mh2),t1(t-_mh2),u1(u-_mh2);
  return mf2*4.*sqrt(2.*s*t*u)*(-4.*(1./(u*t)+1./(u*u1)+1./(t*t1))
				-4.*((2.*s+t)*_bi[k]/sqr(u1)+(2.*s+u)*_bi[j]/sqr(t1))/s
				-(s-4.*mf2)*(s1*_ci[i1]+(u-s)*_ci[j1]+(t-s)*_ci[k1])/(s*t*u)
				-8.*mf2*(_ci[j1]/(t*t1)+_ci[k1]/(u*u1))
				+0.5*(s-4.*mf2)*(s*t*_di[k]+u*s*_di[j]-u*t*_di[i])/(s*t*u)
				+4.*mf2*_di[i]/s
				-2.*(u*_ci[k]+t*_ci[j]+u1*_ci[k1]+t1*_ci[j1]-u*t*_di[i])/sqr(s));
}

complex<Energy> GGtoHHardGenerator::me2(Energy2 s,Energy2 t,Energy2 u,
				       Energy2 mf2) const {
  Energy2 s1(s-_mh2),t1(t-_mh2),u1(u-_mh2);
  return mf2*4.*sqrt(2.*s*t*u)*(4.*_mh2+(_mh2-4.*mf2)*(s1*_ci[4]+t1*_ci[5]+u1*_ci[6])
				-0.5*(_mh2-4.*mf2)*(s*t*_di[3]+u*s*_di[2]+u*t*_di[1]) )/
    (s*t*u);
}

Complex GGtoHHardGenerator::F(double x) {
  if(x<.25) {
    double root = sqrt(1.-4.*x);
    Complex pii(0.,Constants::pi);
    return 0.5*sqr(log((1.+root)/(1.-root))-pii);
  }
  else {
    return -2.*sqr(asin(0.5/sqrt(x)));
  }
}

Energy4 GGtoHHardGenerator::loME() {
  Complex I(0);
  if(_massopt==0) {
    for ( int ix=_minloop; ix<=_maxloop; ++ix ) {
      double x = sqr(getParticleData(ix)->mass())/_mh2;
      I += 3.*x*(2.+(4.*x-1.)*F(x));
    }
  }
  else {
    I = 1.;
  }
  return sqr(_mh2)/576./Constants::pi*norm(I);
}

tPDPtr GGtoHHardGenerator::quarkFlavour(tcPDFPtr pdf, Energy2 scale, 
				       double x, tcBeamPtr beam,
				       double & pdfweight, bool anti) {
  vector<double> weights;
  vector<tPDPtr> partons;
  pdfweight = 0.;
  if(!anti) {
    for(int ix=1;ix<=5;++ix) {
      partons.push_back(getParticleData(ix));
      weights.push_back(pdf->xfx(beam,partons.back(),scale,x));
      pdfweight += weights.back();
    }
  }
  else {
    for(int ix=1;ix<=5;++ix) {
      partons.push_back(getParticleData(-ix));
      weights.push_back(pdf->xfx(beam,partons.back(),scale,x));
      pdfweight += weights.back();
    }
  }
  double wgt=UseRandom::rnd()*pdfweight;
  for(unsigned int ix=0;ix<weights.size();++ix) {
    if(wgt<=weights[ix]) return partons[ix];
    wgt -= weights[ix];
  }
  assert(false);
  return tPDPtr();
}


Energy2 GGtoHHardGenerator::ggME(Energy2 s, Energy2 t, Energy2 u) {
  Energy2 output;
  if(_massopt==0) {
    complex<Energy> me[2][2][2];
    me[1][1][1] = ZERO;
    me[1][1][0] = ZERO;
    me[0][1][0] = ZERO;
    me[0][1][1] = ZERO;
    for(int ix=_minloop;ix<=_maxloop;++ix) {
      Energy2 mf2=sqr(getParticleData(ix)->mass());
      _bi[1]=B(s,mf2);
      _bi[2]=B(u,mf2);
      _bi[3]=B(t,mf2);
      _bi[4]=B(_mh2,mf2);
      _bi[1]=_bi[1]-_bi[4];
      _bi[2]=_bi[2]-_bi[4];
      _bi[3]=_bi[3]-_bi[4];
      _ci[1]=C(s,mf2);
      _ci[2]=C(u,mf2);
      _ci[3]=C(t,mf2);
      _ci[7]=C(_mh2,mf2);
      _ci[4]=(s*_ci[1]-_mh2*_ci[7])/(s-_mh2);
      _ci[5]=(u*_ci[2]-_mh2*_ci[7])/(u-_mh2);
      _ci[6]=(t*_ci[3]-_mh2*_ci[7])/(t-_mh2);
      _di[1]=D(t,u,s,mf2);
      _di[2]=D(s,t,u,mf2);
      _di[3]=D(s,u,t,mf2);
      me[1][1][1]+=me1(s,u,t,mf2,1,2,3,4,5,6);
      me[1][1][0]+=me2(s,u,t,mf2);
      me[0][1][0]+=me1(u,s,t,mf2,2,1,3,5,4,6);
      me[0][1][1]+=me1(t,u,s,mf2,3,2,1,6,5,4);
    }
    me[0][0][0]=-me[1][1][1];
    me[0][0][1]=-me[1][1][0];
    me[1][0][1]=-me[0][1][0];
    me[1][0][0]=-me[0][1][1];
    output = real(me[0][0][0]*conj(me[0][0][0])+
		  me[0][0][1]*conj(me[0][0][1])+
		  me[0][1][0]*conj(me[0][1][0])+
		  me[0][1][1]*conj(me[0][1][1])+
		  me[1][0][0]*conj(me[1][0][0])+
		  me[1][0][1]*conj(me[1][0][1])+
		  me[1][1][0]*conj(me[1][1][0])+
		  me[1][1][1]*conj(me[1][1][1]));
    output *= 3./8.;
  }
  else {
    output=32./3.*
      (pow<4,1>(s)+pow<4,1>(t)+pow<4,1>(u)+pow<4,1>(_mh2))/s/t/u;
  }
  // spin and colour factors
  return output/4./64.;
}

Energy2 GGtoHHardGenerator::qgME(Energy2 s, Energy2 t, Energy2 u) {
  Energy2 output;
  if(_massopt==0) {
    complex<Energy2> A(ZERO);
    Energy2 si(u-_mh2);
    for(int ix=_minloop;ix<=_maxloop;++ix) {
      Energy2 mf2=sqr(getParticleData(ix)->mass());
      A += mf2*(2.+2.*double(u/si)*(B(u,mf2)-B(_mh2,mf2))
 		+double((4.*mf2-s-t)/si)*Complex(u*C(u,mf2)-_mh2*C(_mh2,mf2)));
    }
    output =-4.*(sqr(s)+sqr(t))/sqr(si)/u*real(A*conj(A));
  }
  else{
    output =-4.*(sqr(s)+sqr(t))/u/9.;
  }
  // final colour/spin factors
  return output/24.;
}

Energy2 GGtoHHardGenerator::qbargME(Energy2 s, Energy2 t, Energy2 u) {
  Energy2 output;
  if(_massopt==0) {
    complex<Energy2> A(ZERO);
    Energy2 si(u-_mh2);
    for(int ix=_minloop;ix<=_maxloop;++ix) {
      Energy2 mf2=sqr(getParticleData(ix)->mass());
      A+=mf2*(2.+2.*double(u/si)*(B(u,mf2)-B(_mh2,mf2))
	      +double((4.*mf2-s-t)/si)*Complex(u*C(u,mf2)-_mh2*C(_mh2,mf2)));
    }
    output =-4.*(sqr(s)+sqr(t))/sqr(si)/u*real(A*conj(A));
  }
  else {
    output =-4.*(sqr(s)+sqr(t))/u/9.;
  }
  // final colour/spin factors
  return output/24.;
}

double GGtoHHardGenerator::getResult(int emis_type, Energy pt, double yj,
				     tcPDPtr & outParton) {
  Energy2 s=sqr(generator()->maximumCMEnergy());
  Energy2 scale = _mh2+sqr(pt);
  Energy  et=sqrt(scale);  // Think before you go moving this!
  scale *= sqr(scaleFact_);
  // longitudinal real correction fractions
  double x  = pt*exp( yj)/sqrt(s)+et*exp( _yh)/sqrt(s);
  double y  = pt*exp(-yj)/sqrt(s)+et*exp(-_yh)/sqrt(s);
  // reject if outside region
  if(x<0.||x>1.||y<0.||y>1.||x*y<_mh2/s) return 0.;
  // longitudinal born fractions
  double x1 = _mass*exp( _yh)/sqrt(s);          
  double y1 = _mass*exp(-_yh)/sqrt(s);
  // mandelstam variables
  Energy2 th = -sqrt(s)*x*pt*exp(-yj);
  Energy2 uh = -sqrt(s)*y*pt*exp( yj);
  Energy2 sh = _mh2-th-uh;
  InvEnergy2 res = InvEnergy2();
  // pdf part of the cross section
  double pdf[4];
  pdf[0]=_beams[0]->pdf()->xfx(_beams[0],_partons[0],_mh2,x1);
  pdf[1]=_beams[1]->pdf()->xfx(_beams[1],_partons[1],_mh2,y1);
  // g g -> H g
  if(emis_type==0) {
    outParton = _partons[1];
    pdf[2]=_beams[0]->pdf()->xfx(_beams[0],_partons[0],scale,x);
    pdf[3]=_beams[1]->pdf()->xfx(_beams[1],_partons[1],scale,y);
    res = ggME(sh,uh,th)/loME();
  }
  // q g -> H q 
  else if(emis_type==1) {
    outParton = quarkFlavour(_beams[0]->pdf(),scale,x,_beams[0],pdf[2],false);
    pdf[3]=_beams[1]->pdf()->xfx(_beams[1],_partons[1],scale,y);
    res = qgME(sh,uh,th)/loME();
  }
  // g q -> H q
  else if(emis_type==2) {
    pdf[2]=_beams[0]->pdf()->xfx(_beams[0],_partons[0],scale,x);
    outParton = quarkFlavour(_beams[1]->pdf(),scale,y,_beams[1],pdf[3],false);
    res = qgME(sh,th,uh)/loME();
  }
  // qbar g -> H qbar
  else if(emis_type==3) {
    outParton = quarkFlavour(_beams[0]->pdf(),scale,x,_beams[0],pdf[2],true);
    pdf[3]=_beams[1]->pdf()->xfx(_beams[1],_partons[1],scale,y);
    res = qbargME(sh,uh,th)/loME();
  }
  // g qbar -> H qbar
  else if(emis_type==4) {
    pdf[2]=_beams[0]->pdf()->xfx(_beams[0],_partons[0],scale,x);
    outParton = quarkFlavour(_beams[1]->pdf(),scale,y,_beams[1],pdf[3],true);
    res = qbargME(sh,th,uh)/loME();
  }
  //deals with pdf zero issue at large x
  if(pdf[0]<=0.||pdf[1]<=0.||pdf[2]<=0.||pdf[3]<=0.) {
    res = ZERO;
  }
  else {
    res *= pdf[2]*pdf[3]/pdf[0]/pdf[1]*_mh2/sh;
  }
  double the_result(0.);
  if(alphaScale_==0) {
    the_result=_alphaS->ratio(sqr(pt*scaleFact_))
      /8./sqr(Constants::pi)*_mh2/sh*GeV*pt*res;
  } else {
    the_result=_alphaS->ratio(scale)
      /8./sqr(Constants::pi)*_mh2/sh*GeV*pt*res;
  }
  return the_result;
}
