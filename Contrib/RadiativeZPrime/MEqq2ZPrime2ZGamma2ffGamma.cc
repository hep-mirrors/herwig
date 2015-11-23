// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEqq2ZPrime2ZGamma2ffGamma class.
//

#include "MEqq2ZPrime2ZGamma2ffGamma.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "RadiativeZPrimeModel.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace RadiativeZPrime;

void MEqq2ZPrime2ZGamma2ffGamma::doinit() {
  MEBase::doinit(); 
  _z0     = getParticleData(ParticleID::Z0);
  _zPrime = getParticleData(32); 
  tcSMPtr sm = generator()->standardModel();
  tcRadiativeZPrimeModelPtr model = 
    dynamic_ptr_cast<tcRadiativeZPrimeModelPtr>(generator()->standardModel());
  if(!model) throw Exception() << "Must be using the RadiativeZPrimeModel in "
			       << "MEqq2ZPrime2ZGamma::doinit()" << Exception::abortnow;
  _theFFZPrimeVertex     = model->vertexFFZPrime();
  _theGammaZPrimeZVertex = model->vertexGammaZPrimeZ();
  _theFFZVertex          = model->vertexFFZ();
}

void MEqq2ZPrime2ZGamma2ffGamma::getDiagrams() const {
  // find possible Z decays
  typedef Selector<tDMPtr> DecaySelector;
  DecaySelector Zdec = _z0->decaySelector();
  vector<PDPair> Zdecays;
  for(DecaySelector::const_iterator cit=Zdec.begin();cit!=Zdec.end();++cit) {
    if(cit->second->orderedProducts().size()!=2) continue;
    if(cit->second->orderedProducts()[0]->id()>0)
      Zdecays.push_back(make_pair(cit->second->orderedProducts()[0],
				  cit->second->orderedProducts()[1]));
    else
      Zdecays.push_back(make_pair(cit->second->orderedProducts()[1],
				  cit->second->orderedProducts()[0]));
  }
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  for(unsigned int i = 1; i <= _maxflavour; ++i) {
    tcPDPtr q  = getParticleData(long(i));
    tcPDPtr qb = q->CC();
    for(unsigned int iz=0;iz<Zdecays.size();++iz) {
      add(new_ptr((Tree2toNDiagram(2), q, qb, 1, _zPrime, 3, gamma, 3, _z0, 
		   5, Zdecays[iz].first,5, Zdecays[iz].second, -1)));
    }
  }
}

Energy2 MEqq2ZPrime2ZGamma2ffGamma::scale() const {
  return sHat();
}

int MEqq2ZPrime2ZGamma2ffGamma::nDim() const {
  return 5;
}

bool MEqq2ZPrime2ZGamma2ffGamma::generateKinematics(const double * r) {
  // initialize jacobian
  jacobian(1.);
  // cms energy
  Energy ecm=sqrt(sHat());
  // first generate the mass of the off-shell gauge boson
  // minimum mass of the 
  tcPDVector ptemp;
  ptemp.push_back(mePartonData()[3]);
  ptemp.push_back(mePartonData()[4]);
  Energy2 minMass2 = max(lastCuts().minSij(mePartonData()[3],mePartonData()[4]),
			 lastCuts().minS(ptemp));
  // minimum pt of the jet
  Energy ptmin = max(lastCuts().minKT(mePartonData()[2]),
		     lastCuts().minKT(_z0));
  // maximum mass of the gauge boson so pt is possible
  Energy2 maxMass2 = min(ecm*(ecm-2.*ptmin),lastCuts().maxS(ptemp));
  if(maxMass2<=ZERO||minMass2<ZERO) return false;
  // also impose the limits from the ParticleData object
  minMass2 = max(minMass2,sqr(_z0->massMin()));
  maxMass2 = min(maxMass2,sqr(_z0->massMax()));
  // also impose the limits from the ParticleData object
  if(maxMass2<minMass2) return false;
  // generation of the mass
  Energy  M(_z0->mass()),Gamma(_z0->width());
  Energy2 M2(sqr(M)),MG(M*Gamma);
  double rhomin = atan2((minMass2-M2),MG);
  double rhomax = atan2((maxMass2-M2),MG);
  _mz2=M2+MG*tan(rhomin+r[1]*(rhomax-rhomin));
  Energy mz=sqrt(_mz2);
  InvEnergy2 emjac = MG/(rhomax-rhomin)/(sqr(_mz2-M2)+sqr(MG));
  // jacobian
  jacobian(jacobian()/sHat()/emjac);
  // set the masses of the outgoing particles to 2-2 scattering
  meMomenta()[2].setMass(ZERO);
  Lorentz5Momentum pz(mz);
  // generate the polar angle of the hard scattering
  double ctmin(-1.0), ctmax(1.0);
  Energy q(ZERO);
  try {
    q = SimplePhaseSpace::getMagnitude(sHat(), meMomenta()[2].mass(),mz);
  } 
  catch ( ImpossibleKinematics ) {
    return false;
  }	    
  Energy2 pq = sqrt(sHat())*q;
  if ( ptmin > ZERO ) {
    double ctm = 1.0 - sqr(ptmin/q);
    if ( ctm <= 0.0 ) return false;
    ctmin = max(ctmin, -sqrt(ctm));
    ctmax = min(ctmax,  sqrt(ctm));
  }
  if ( ctmin >= ctmax ) return false;
  double cth = getCosTheta(ctmin, ctmax, r[0]);
  Energy pt  = q*sqrt(1.0-sqr(cth));
  double phi = 2.0*Constants::pi*r[2];
  meMomenta()[2].setVect(Momentum3( pt*sin(phi), pt*cos(phi), q*cth));
  pz.setVect(            Momentum3(-pt*sin(phi),-pt*cos(phi),-q*cth));
  meMomenta()[2].rescaleEnergy();
  pz.rescaleEnergy();
  // generate the momenta of the Z decay products
  meMomenta()[3].setMass(mePartonData()[3]->mass());
  meMomenta()[4].setMass(mePartonData()[4]->mass());
  Energy q2 = ZERO;
  try {
    q2 = SimplePhaseSpace::getMagnitude(_mz2, meMomenta()[3].mass(),
					meMomenta()[4].mass());
  } catch ( ImpossibleKinematics ) {
    return false;
  }
  double cth2 =-1.+2.*r[3];
  double phi2=Constants::twopi*r[4];
  Energy pt2 =q2*sqrt(1.-sqr(cth2));
  Lorentz5Momentum pl[2]={Lorentz5Momentum( pt2*cos(phi2), pt2*sin(phi2), q2*cth2,ZERO,
					    meMomenta()[3].mass()),
			  Lorentz5Momentum(-pt2*cos(phi2),-pt2*sin(phi2),-q2*cth2,ZERO,
					   meMomenta()[4].mass())};
  pl[0].rescaleEnergy();
  pl[1].rescaleEnergy();
  Boost boostv(pz.boostVector());
  pl[0].boost(boostv);
  pl[1].boost(boostv);
  meMomenta()[3] = pl[0];
  meMomenta()[4] = pl[1];
  // check passes all the cuts
  vector<LorentzMomentum> out(3);
  tcPDVector tout(3);
  for(unsigned int ix=0;ix<3;++ix) {
    out[ ix] = meMomenta()[ix+2];
    tout[ix] = mePartonData()[ix+2];
  }
  if ( !lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]) )
    return false;
  // jacobian
  jacobian((pq/sHat())*Constants::pi*jacobian()/8./sqr(Constants::pi)*q2/mz);
  return true;
}

double MEqq2ZPrime2ZGamma2ffGamma::me2() const {
  InvEnergy2 output(ZERO);
  // construct spinors for the leptons
  vector<SpinorBarWaveFunction> lm;
  vector<SpinorWaveFunction>    lp;
  SpinorBarWaveFunction lmout(meMomenta()[3],mePartonData()[3],outgoing);
  SpinorWaveFunction    lpout(meMomenta()[4],mePartonData()[4],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    lmout.reset(ix);lm.push_back(lmout);
    lpout.reset(ix);lp.push_back(lpout);
  }
  vector<SpinorWaveFunction>     fin;
  vector<SpinorBarWaveFunction>  ain;
  vector<VectorWaveFunction> gout;
  SpinorWaveFunction    qin (meMomenta()[0],mePartonData()[0],incoming);
  SpinorBarWaveFunction qbin(meMomenta()[1],mePartonData()[1],incoming);
  VectorWaveFunction   glout(meMomenta()[2],mePartonData()[2],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    qin.reset(ix)  ;  fin.push_back(qin);
    qbin.reset(ix) ;  ain.push_back(qbin);
    glout.reset(2*ix); gout.push_back(glout);
  }
  output=qqbarME(fin,ain,gout,lm,lp);
  return output*sHat();
}

CrossSection MEqq2ZPrime2ZGamma2ffGamma::dSigHatDR() const {
  return me2()*jacobian()/(16.0*sqr(Constants::pi)*sHat())*sqr(hbarc);
}

unsigned int MEqq2ZPrime2ZGamma2ffGamma::orderInAlphaS() const {
  return 0;
}

unsigned int MEqq2ZPrime2ZGamma2ffGamma::orderInAlphaEW() const {
  return 4;
}

Selector<MEBase::DiagramIndex>
MEqq2ZPrime2ZGamma2ffGamma::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) sel.insert(1., i);
  return sel;
}

Selector<const ColourLines *>
MEqq2ZPrime2ZGamma2ffGamma::colourGeometries(tcDiagPtr diag) const {
  static const ColourLines c1=ColourLines("1 -2");
  static const ColourLines c2=ColourLines("1 -2,  6 -7");
  Selector<const ColourLines *> sel;
  if(mePartonData()[3]->coloured()) sel.insert(1.0,&c2);
  else                              sel.insert(1.0,&c1);
  return sel;
}

void MEqq2ZPrime2ZGamma2ffGamma::persistentOutput(PersistentOStream & os) const {
  os << _theFFZVertex << _theFFZPrimeVertex << _theGammaZPrimeZVertex
     << _z0 << _zPrime << _maxflavour;
}

void MEqq2ZPrime2ZGamma2ffGamma::persistentInput(PersistentIStream & is, int) {
  is >> _theFFZVertex >> _theFFZPrimeVertex >> _theGammaZPrimeZVertex
     >> _z0 >> _zPrime >> _maxflavour;
}

ClassDescription<MEqq2ZPrime2ZGamma2ffGamma> MEqq2ZPrime2ZGamma2ffGamma::initMEqq2ZPrime2ZGamma2ffGamma;
// Definition of the static class description member.

void MEqq2ZPrime2ZGamma2ffGamma::Init() {

  static ClassDocumentation<MEqq2ZPrime2ZGamma2ffGamma> documentation
    ("There is no documentation for the MEqq2ZPrime2ZGamma2ffGamma class");

  static Parameter<MEqq2ZPrime2ZGamma2ffGamma,unsigned int> interfaceMaxFlavour
    ("MaxFlavour",
     "The heaviest incoming quark flavour this matrix element is allowed to handle",
     &MEqq2ZPrime2ZGamma2ffGamma::_maxflavour, 5, 1, 6,
     false, false, Interface::limited);

}

double MEqq2ZPrime2ZGamma2ffGamma::
getCosTheta(double ctmin, double ctmax, const double r) {
  double cth = 0.0;
  double zmin = 0.5*(1.0 - ctmax);
  double zmax = 0.5*(1.0 - ctmin);
  if ( zmin <= 0.0 || zmax >= 1.0 ) {
    jacobian((ctmax - ctmin)*jacobian());
    cth = ctmin + r*(ctmax-ctmin);
  } else {
    double A1 = (2.0*zmax - 1.0)/(zmax*(1.0-zmax));
    double A0 = (2.0*zmin - 1.0)/(zmin*(1.0-zmin));
    double A = r*(A1 - A0) + A0;
    double z = A < 2.0? 2.0/(sqrt(sqr(A) + 4.0) + 2 - A):
      0.5*(A - 2.0 + sqrt(sqr(A) + 4.0))/A;
    cth = 1.0 - 2.0*z;
    jacobian((2.0*(A1 - A0)*sqr(z)*sqr(1.0 - z)/(sqr(z) + sqr(1.0 - z)))*jacobian());
  }
  return cth;
}  

InvEnergy2 MEqq2ZPrime2ZGamma2ffGamma::qqbarME(vector<SpinorWaveFunction> & fin,
					       vector<SpinorBarWaveFunction> & ain,
					       vector<VectorWaveFunction> & gout,
					       vector<SpinorBarWaveFunction> & lm,
					       vector<SpinorWaveFunction> & lp,
					       bool calc) const {
  // scale
  Energy2 mb2(scale());
  // if calculation spin corrections construct the me
  if(calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
					     PDT::Spin1,PDT::Spin1Half,
					     PDT::Spin1Half));
  // some integers
  unsigned int ihel1,ihel2,ohel1,ohel2,ohel3;
  // compute the leptonic Z currents for speed
  VectorWaveFunction bcurr[2][2];
  for(ohel2=0;ohel2<2;++ohel2) {
    for(ohel3=0;ohel3<2;++ohel3) {
      bcurr[ohel2][ohel3]=
	_theFFZVertex->evaluate(_mz2,1,_z0,lp[ohel3],lm[ohel2]);
    }
  }

  // compute the matrix elements
  double me=0.;
  Complex diag;
  VectorWaveFunction inter;
  for(ihel1=0;ihel1<2;++ihel1) {
    for(ihel2=0;ihel2<2;++ihel2) {
      // intermediate for Z'
      inter=_theFFZPrimeVertex->evaluate(mb2,1,_zPrime,fin[ihel1],ain[ihel2]);
      for(ohel1=0;ohel1<2;++ohel1) {
	for(ohel2=0;ohel2<2;++ohel2) {
	  for(ohel3=0;ohel3<2;++ohel3) {
	    diag = _theGammaZPrimeZVertex->evaluate(mb2,gout[ohel1],inter,
						    bcurr[ohel2][ohel3]);
	    me += norm(diag);
	    if(calc) _me(ihel1,ihel2,2*ohel1,ohel2,ohel3) = diag;
	  }
	}
      }
    }
  }
  // results
  // initial state spin and colour average
  double colspin = 1./3./4.;
  // and for Z decay products
  if(mePartonData()[3]->coloured()) colspin *= 3.;
  me *= colspin;
  return me*UnitRemoval::InvE2;
}

void  MEqq2ZPrime2ZGamma2ffGamma::constructVertex(tSubProPtr sub) {
//   // extract the particles in the hard process
//   ParticleVector hard(5);
//   // incoming
//   hard[0]=sub->incoming().first;
//   hard[1]=sub->incoming().second;
//   if(hard[0]->id()<0) swap(hard[0],hard[1]);
//   // outgoing
//   for(unsigned int ix=0;ix<3;++ix) {
//     unsigned int iloc;
//     PPtr mother=sub->outgoing()[ix]->parents()[0];
//     if(mother&&(mother->id()==ParticleID::gamma||mother->id()==ParticleID::Z0)) {
//       if(sub->outgoing()[ix]->id()>0) iloc=3;
//       else                            iloc=4;
//     }
//     else iloc=2;
//     hard[iloc]=sub->outgoing()[ix];
//   }
//   // wavefunctions for the Z decay products
//   vector<SpinorBarWaveFunction> lm;
//   vector<SpinorWaveFunction>    lp;
//   SpinorBarWaveFunction(lm,hard[3],outgoing,true,true);
//   SpinorWaveFunction   (lp,hard[4],outgoing,true,true);
//   vector<SpinorWaveFunction>     fin;
//   vector<SpinorBarWaveFunction>  ain;
//   vector<VectorWaveFunction> gout;
//   SpinorWaveFunction   (fin ,hard[0],incoming,false,true);
//   SpinorBarWaveFunction(ain ,hard[1],incoming,false,true);
//   VectorWaveFunction   (gout,hard[2],outgoing,true ,true,true);
//   gout[1]=gout[2];
//   qqbarME(fin,ain,gout,lm,lp,true);
//   // construct the vertex
//   HardVertexPtr hardvertex=new_ptr(HardVertex());
//   // set the matrix element for the vertex
//   hardvertex->ME(_me);
//   // set the pointers and to and from the vertex
//   for(unsigned int ix=0;ix<5;++ix)
//     (hard[ix]->spinInfo())->productionVertex(hardvertex);
}
