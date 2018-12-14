// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEfftoVH class.
//

#include "MEfftoVH.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/PDF/PolarizedBeamParticleData.h"

using namespace Herwig;

void MEfftoVH::persistentOutput(PersistentOStream & os) const {
  os << _shapeopt << _wplus << _wminus << _z0 << _higgs
     << _vertexFFW << _vertexFFZ << _vertexWWH << _maxflavour
     << ounit(_mh,GeV) << ounit(_wh,GeV) << _hmass;
}

void MEfftoVH::persistentInput(PersistentIStream & is, int) {
  is >> _shapeopt >> _wplus >> _wminus >> _z0 >> _higgs 
     >> _vertexFFW >> _vertexFFZ >> _vertexWWH >> _maxflavour
     >> iunit(_mh,GeV) >> iunit(_wh,GeV) >> _hmass;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<MEfftoVH,DrellYanBase>
describeHerwigMEfftoVH("Herwig::MEfftoVH", "Herwig.so");

void MEfftoVH::Init() {

  static ClassDocumentation<MEfftoVH> documentation
    ("The MEfftoVH class is the base class for the Bjirken process f fbar -> V H");

  static Switch<MEfftoVH,unsigned int> interfaceShapeOption
    ("ShapeScheme",
     "Option for the treatment of the Higgs resonance shape",
     &MEfftoVH::_shapeopt, 2, false, false);
  static SwitchOption interfaceStandardShapeFixed
    (interfaceShapeOption,
     "FixedBreitWigner",
     "Breit-Wigner s-channel resonanse",
     1);
  static SwitchOption interfaceStandardShapeRunning
    (interfaceShapeOption,
     "MassGenerator",
     "Use the mass generator to give the shape",
     2);
  static SwitchOption interfaceStandardShapeOnShell
    (interfaceShapeOption,
     "OnShell",
     "Produce an on-shell Higgs boson",
     0);

  static Parameter<MEfftoVH,unsigned int> interfaceMaxFlavour
    ( "MaxFlavour",
      "The heaviest incoming quark flavour this matrix element is allowed to handle "
      "(if applicable).",
      &MEfftoVH::_maxflavour, 5, 1, 5, false, false, true);

}

unsigned int MEfftoVH::orderInAlphaS() const {
  return 0;
}

unsigned int MEfftoVH::orderInAlphaEW() const {
  return 3;
}

Energy2 MEfftoVH::scale() const {
  return sHat();
}

int MEfftoVH::nDim() const {
  return 4 + int(_shapeopt>0);
}

void MEfftoVH::setKinematics() {
  DrellYanBase::setKinematics();
}

Selector<MEBase::DiagramIndex>
MEfftoVH::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    sel.insert(1.0, i);
  return sel;
}

Selector<const ColourLines *>
MEfftoVH::colourGeometries(tcDiagPtr ) const {
  static ColourLines c1("");
  static ColourLines c2("6 -7");
  static ColourLines c3("1 -2");
  static ColourLines c4("1 -2, 6 -7");
  Selector<const ColourLines *> sel;
  if(mePartonData()[0]->coloured()) {
    if(mePartonData()[4]->coloured()) sel.insert(1.0, &c4);
    else                              sel.insert(1.0, &c3);
  }
  else {
    if(mePartonData()[4]->coloured()) sel.insert(1.0, &c2);
    else                              sel.insert(1.0, &c1);
  }
  return sel;
}

void MEfftoVH::doinit() {
  DrellYanBase::doinit();
  // get the vedrtex pointers from the SM object
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(hwsm) {
    _vertexFFW = hwsm->vertexFFW();
    _vertexFFZ = hwsm->vertexFFZ();
  }
  else throw InitException() << "Wrong type of StandardModel object in "
			     << "MEfftoVH::doinit() the Herwig"
			     << " version must be used" 
			     << Exception::runerror;
  // get the particle data objects for the intermediates
  _wplus  = getParticleData(ParticleID::Wplus );
  _wminus = getParticleData(ParticleID::Wminus);
  _z0     = getParticleData(ParticleID::Z0);
  // higgs stuff
  _mh = _higgs->mass();
  _wh = _higgs->width();
  if(_higgs->massGenerator()) {
    _hmass=dynamic_ptr_cast<GenericMassGeneratorPtr>(_higgs->massGenerator());
  }
  if(_shapeopt==2&&!_hmass) 
    throw InitException()
      << "If using the mass generator for the line shape in MEfftoVH::doinit()"
      << "the mass generator must be an instance of the GenericMassGenerator class"
      << Exception::runerror;
}

double MEfftoVH::me2() const {
  vector<SpinorWaveFunction>    fin,aout;
  vector<SpinorBarWaveFunction> ain,fout;
  SpinorWaveFunction       q(meMomenta()[0],mePartonData()[0],incoming);
  SpinorBarWaveFunction qbar(meMomenta()[1],mePartonData()[1],incoming);
  SpinorBarWaveFunction    f(meMomenta()[3],mePartonData()[3],outgoing);
  SpinorWaveFunction    fbar(meMomenta()[4],mePartonData()[4],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    q.reset(ix)   ; fin.push_back(q);
    qbar.reset(ix); ain.push_back(qbar);
    f.reset(ix)   ;fout.push_back(f);
    fbar.reset(ix);aout.push_back(fbar);
  }
  return helicityME(fin,ain,fout,aout,false)*sHat()*UnitRemoval::InvE2;
}

double MEfftoVH::helicityME(vector<SpinorWaveFunction>    & fin ,
			    vector<SpinorBarWaveFunction> & ain ,
			    vector<SpinorBarWaveFunction> & fout,
			    vector<SpinorWaveFunction>    & aout,
			    bool calc) const {
  // scale
  Energy2 mb2(scale());
  // matrix element to be stored
  ProductionMatrixElement menew(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0,
				PDT::Spin1Half,PDT::Spin1Half);
  // work out the id of the vector boson 
  int incharge = mePartonData()[0]->iCharge()+mePartonData()[1]->iCharge();
  tcPDPtr vec;
  if(incharge==0)     vec = _z0;
  else if(incharge>0) vec = _wplus;
  else                vec = _wminus;
  // vertex for vector boson
  AbstractFFVVertexPtr vertex = vec==_z0 ? _vertexFFZ : _vertexFFW;
  // wavefunction for the scalar
  ScalarWaveFunction higgs(meMomenta()[2],mePartonData()[2],1.,outgoing);
  // calculate the matrix element
  VectorWaveFunction inter[2];
  unsigned int ihel1,ihel2,ohel1,ohel2;
  Complex diag;
  double me(0.);
  for(ihel1=0;ihel1<2;++ihel1) {
    for(ihel2=0;ihel2<2;++ihel2) {
      // wavefunction for the intermediate 1st vector
      inter[0] = vertex->evaluate(mb2,1,vec,fin[ihel1],ain[ihel2]);
      // boson decay piece
      for(ohel1=0;ohel1<2;++ohel1) {
	for(ohel2=0;ohel2<2;++ohel2) {
	  inter[1] = vertex->evaluate(sqr(vec->mass()),1,vec,
				      aout[ohel2],fout[ohel1]);
      	  diag = _vertexWWH->evaluate(mb2,inter[1],inter[0],higgs);
      	  me += norm(diag);
	  menew(ihel1,ihel2,0,ohel1,ohel2) = diag;
	}
      }
    }
  }
  // spin factor
  me *=0.25;
  tcPolarizedBeamPDPtr beam[2] = 
    {dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[0]),
     dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[1])};
  if( beam[0] || beam[1] ) {
    RhoDMatrix rho[2] = 
      {beam[0] ? beam[0]->rhoMatrix() : RhoDMatrix(mePartonData()[0]->iSpin()),
       beam[1] ? beam[1]->rhoMatrix() : RhoDMatrix(mePartonData()[1]->iSpin())};
    me = menew.average(rho[0],rho[1]);
  }
  // incoming colour factor
  if(mePartonData()[0]->coloured()) me /= 3.;
  // outgoing colour factor
  if(mePartonData()[3]->coloured()) me *= 3.;
  if(calc) _me.reset(menew);
  return me;
}

void MEfftoVH::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
  hard.push_back(sub->outgoing()[2]);
  // ensure right order
  if(hard[0]->id()<0) swap(hard[0],hard[1]);
  if(hard[3]->dataPtr()->iSpin()==PDT::Spin0) swap(hard[2],hard[3]);
  if(hard[4]->dataPtr()->iSpin()==PDT::Spin0) swap(hard[2],hard[4]);
  if(hard[3]->id()<0) swap(hard[3],hard[4]);
  vector<SpinorWaveFunction>    fin,aout;
  vector<SpinorBarWaveFunction> ain,fout;
  SpinorWaveFunction(   fin ,hard[0],incoming,false,true);
  SpinorBarWaveFunction(ain ,hard[1],incoming,false,true);
  ScalarWaveFunction(        hard[2],outgoing,true);
  SpinorBarWaveFunction(fout,hard[3],outgoing,true ,true);
  SpinorWaveFunction(   aout,hard[4],outgoing,true ,true);
  helicityME(fin,ain,fout,aout,true);
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<5;++ix) {
    tcSpinPtr spin = hard[ix]->spinInfo();
    if(ix<2) {
      tcPolarizedBeamPDPtr beam = 
	dynamic_ptr_cast<tcPolarizedBeamPDPtr>(hard[ix]->dataPtr());
      if(beam) spin->rhoMatrix() = beam->rhoMatrix();
    }
    spin->productionVertex(hardvertex);
  }
}

bool MEfftoVH::generateKinematics(const double * r) {
  using Constants::pi;
  // workout the ID of the vector boson
  tcPDPtr vec = mePartonData()[0]->iCharge()+mePartonData()[1]->iCharge()!=0
    ? _wplus : _z0;
  // order determined randomly
  Energy e = sqrt(sHat())/2.0;
  Energy mh(_mh),mv;
  double jac(1.);
  if(UseRandom::rndbool()) {
    double rhomax,rhomin;
    // generate the mass of the Higgs
    if(_shapeopt!=0) {
      Energy mhmax = min(2.*e-vec->massMin(),mePartonData()[2]->massMax());
      Energy mhmin = max(ZERO               ,mePartonData()[2]->massMin());
      if(mhmax<=mhmin) return false;
      rhomin = atan2((sqr(mhmin)-sqr(_mh)), _mh*_wh);
      rhomax = atan2((sqr(mhmax)-sqr(_mh)), _mh*_wh);
      mh = sqrt(_mh*_wh*tan(rhomin+r[3]*(rhomax-rhomin))+sqr(_mh));
      jac *= rhomax-rhomin;
    }
    // generate the mass of the vector boson
    Energy2 mvmax2 = sqr(min(2.*e-mh,vec->massMax()));
    Energy2 mvmin2 = sqr(vec->massMin());
    if(mvmax2<=mvmin2) return false; 
    rhomin = atan2((mvmin2-sqr(vec->mass())), vec->mass()*vec->width());
    rhomax = atan2((mvmax2-sqr(vec->mass())), vec->mass()*vec->width());
    mv = sqrt(vec->mass()*vec->width()*tan(rhomin+r[1]*(rhomax-rhomin))
	      +sqr(vec->mass()));
    jac *= rhomax-rhomin;
  }
  else {
    // generate the mass of the vector boson
    Energy2 mvmax2 = sqr(min(2.*e,vec->massMax()));
    Energy2 mvmin2 = sqr(vec->massMin());
    if(mvmax2<=mvmin2) return false; 
    double rhomin = atan2((mvmin2-sqr(vec->mass())), vec->mass()*vec->width());
    double rhomax = atan2((mvmax2-sqr(vec->mass())), vec->mass()*vec->width());
    mv = sqrt(vec->mass()*vec->width()*tan(rhomin+r[1]*(rhomax-rhomin))
	      +sqr(vec->mass()));
    jac *= rhomax-rhomin;
    // generate the mass of the Higgs
    if(_shapeopt!=0) {
      Energy mhmax = min(2.*e-mv,mePartonData()[2]->massMax());
      Energy mhmin = max(ZERO ,mePartonData()[2]->massMin());
      if(mhmax<=mhmin) return false;
      rhomin = atan2((sqr(mhmin)-sqr(_mh)), _mh*_wh);
      rhomax = atan2((sqr(mhmax)-sqr(_mh)), _mh*_wh);
      mh = sqrt(_mh*_wh*tan(rhomin+r[3]*(rhomax-rhomin))+sqr(_mh));
      jac *= rhomax-rhomin;
    }
  }
  if(mh+mv>2.*e) return false;
  // assign masses
  meMomenta()[2].setMass(mh);
  for(unsigned int ix=3;ix<5;++ix) 
    meMomenta()[ix] = Lorentz5Momentum(mePartonData()[ix]->generateMass());
  Energy q = ZERO;
  Lorentz5Momentum pvec(mv);
  try {
    q = SimplePhaseSpace::
      getMagnitude(sHat(), meMomenta()[2].mass(), mv);
  } 
  catch ( ImpossibleKinematics & e) {
    return false;
  }
	
  Energy ptmin = max(lastCuts().minKT(mePartonData()[2]),
   		     lastCuts().minKT(vec));	    
  Energy2 m22  = meMomenta()[2].mass2();
  Energy2 m32  = pvec          .mass2();
  Energy2 e0e2 = 2.0*e*sqrt(sqr(q) + m22);
  Energy2 e1e2 = 2.0*e*sqrt(sqr(q) + m22);
  Energy2 e0e3 = 2.0*e*sqrt(sqr(q) + m32);
  Energy2 e1e3 = 2.0*e*sqrt(sqr(q) + m32);
  Energy2 pq   = 2.0*e*q;
  double ctmin = -1.0,ctmax = 1.0;
  Energy2 thmin = lastCuts().minTij(mePartonData()[0], mePartonData()[2]);
  if ( thmin > ZERO ) ctmax = min(ctmax, (e0e2 - m22 - thmin)/pq);

  thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[2]);
  if ( thmin > ZERO ) ctmin = max(ctmin, (thmin + m22 - e1e2)/pq);

  thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[3]);
  if ( thmin > ZERO ) ctmax = min(ctmax, (e1e3 - m32 - thmin)/pq);

  thmin = lastCuts().minTij(mePartonData()[0], mePartonData()[3]);
  if ( thmin > ZERO ) ctmin = max(ctmin, (thmin + m32 - e0e3)/pq);
  if ( ptmin > ZERO ) {
    double ctm = 1.0 - sqr(ptmin/q);
    if ( ctm <= 0.0 ) return false;
    ctmin = max(ctmin, -sqrt(ctm));
    ctmax = min(ctmax, sqrt(ctm));
  }

  if ( ctmin >= ctmax ) return false;
  jacobian(1.);
  double cth = getCosTheta(ctmin, ctmax, r[0]);
  
  Energy pt = q*sqrt(1.0-sqr(cth));
  double phi = UseRandom::rnd()*2.0*pi;
  meMomenta()[2].setX(pt*sin(phi));
  meMomenta()[2].setY(pt*cos(phi));
  meMomenta()[2].setZ(q*cth);
  
  pvec          .setX(-pt*sin(phi));
  pvec          .setY(-pt*cos(phi));
  pvec          .setZ(-q*cth);
  
  meMomenta()[2].rescaleEnergy();
  pvec          .rescaleEnergy();
  // decay of the vector boson
  bool test=Kinematics::twoBodyDecay(pvec,meMomenta()[3].mass(),
				     meMomenta()[4].mass(),
				     -1.+2*r[2],r[3]*2.*pi,
				     meMomenta()[3],meMomenta()[4]);
  if(!test) return false;
  // check cuts
  vector<LorentzMomentum> out;
  tcPDVector tout;
  for(unsigned int ix=2;ix<5;++ix) {
    out .push_back(meMomenta()[ix]);
    tout.push_back(mePartonData()[ix]);
  }
  if ( !lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]) )
    return false;
  // jacobian factors
  // main piece
  jacobian((pq/sHat())*pi*jacobian());
  // mass piece
  jacobian(jac*jacobian());
  // decay piece
  Energy p2 = Kinematics::pstarTwoBodyDecay(mv,meMomenta()[3].mass(),
					    meMomenta()[4].mass());
  jacobian(p2/mv/8./sqr(pi)*jacobian());
  // jacobian factor for the gauge boson
  jacobian((sqr(sqr(mv)-sqr(vec->mass()))+sqr(vec->mass()*vec->width()))/
	   (vec->mass()*vec->width())*jacobian()/sHat());
  return true;
}

CrossSection MEfftoVH::dSigHatDR() const {
  using Constants::pi;
  // jacobian factor for the higgs
  InvEnergy2 bwfact(ZERO);
  Energy moff =meMomenta()[2].mass();
  if(_shapeopt==1) {
    tcPDPtr h0 = mePartonData()[2]->iSpin()==PDT::Spin0 ?
      mePartonData()[2] : mePartonData()[3];
    bwfact = h0->generateWidth(moff)*moff/pi/
      (sqr(sqr(moff)-sqr(_mh))+sqr(_mh*_wh));
  }
  else if(_shapeopt==2) {
    bwfact = _hmass->BreitWignerWeight(moff);
  }
  double jac1 = _shapeopt!=0 ? 
    double(bwfact*(sqr(sqr(moff)-sqr(_mh))+sqr(_mh*_wh))/(_mh*_wh)) : 1.;
  // answer
  return jac1*me2()*jacobian()/(16.0*sqr(pi)*sHat())*sqr(hbarc);
}
