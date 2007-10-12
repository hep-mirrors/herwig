// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEqq2gZ2ffNason class.
//

#include "MEqq2gZ2ffNason.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "HardVertex.h"

using namespace Herwig;

void MEqq2gZ2ffNason::doinit() throw(InitException) {
  ME2to2Base::doinit();
  _z0=getParticleData(ThePEG::ParticleID::Z0);
  _gamma=getParticleData(ThePEG::ParticleID::gamma);
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(hwsm) {
    _theFFZVertex = hwsm->vertexFFZ();
    _theFFPVertex = hwsm->vertexFFP();
  }
  else
    throw InitException() << "Must be the Herwig++ StandardModel class in "
			  << "MEqq2gZ2ffNason::doinit" << Exception::abortnow;
}

void MEqq2gZ2ffNason::getDiagrams() const {
  // which intermediates to include
  bool gamma = _gammaZ==0 || _gammaZ==1;
  bool Z0    = _gammaZ==0 || _gammaZ==2;
  // loop over the processes we need
  for(unsigned int ix=1;ix<17;++ix) {
    // increment counter to switch between quarks and leptons
    if(ix==7) ix+=4;
    // is it a valid quark process
    bool quark = ix<=6 && (_process==0 || _process==1 || _process-10==ix);
    // is it a valid lepton process
    bool lepton= ix>=11 && ix<=16 
      && (_process==0
	  || _process==2
	  || (_process==3 && ix%2==1)
	  || (_process==4 && ix%2==0)
	  || (ix%2==0 && (ix-10)/2==_process-7)
	  || (ix%2==1 && (ix-9)/2 ==_process-4)
	  );
    // if not a validf process continue
    if(!(quark||lepton)) continue;
    tcPDPtr lm = getParticleData(ix);
    tcPDPtr lp = lm->CC();
    for(unsigned int i = 1; i <= _maxflavour; ++i) {
      tcPDPtr q  = getParticleData(i);
      tcPDPtr qb = q->CC();
      if(Z0)    add(new_ptr((Tree2toNDiagram(2), q, qb, 1, _z0   , 3, lm, 3, lp, -1)));
      if(gamma) add(new_ptr((Tree2toNDiagram(2), q, qb, 1, _gamma, 3, lm, 3, lp, -2)));
    }
  }
}

Energy2 MEqq2gZ2ffNason::scale() const {
  return sHat();
}

double MEqq2gZ2ffNason::me2() const {
  vector<SpinorWaveFunction>    fin,aout;
  vector<SpinorBarWaveFunction> ain,fout;
  SpinorWaveFunction       q(meMomenta()[0],mePartonData()[0],incoming);
  SpinorBarWaveFunction qbar(meMomenta()[1],mePartonData()[1],incoming);
  SpinorBarWaveFunction    f(meMomenta()[2],mePartonData()[2],outgoing);
  SpinorWaveFunction    fbar(meMomenta()[3],mePartonData()[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    q.reset(ix)   ; fin.push_back(q);
    qbar.reset(ix); ain.push_back(qbar);
    f.reset(ix)   ;fout.push_back(f);
    fbar.reset(ix);aout.push_back(fbar);
  }
  return qqbarME(fin,ain,fout,aout,false);
}

unsigned int MEqq2gZ2ffNason::orderInAlphaS() const {
  return 0;
}

unsigned int MEqq2gZ2ffNason::orderInAlphaEW() const {
  return 2;
}

Selector<MEBase::DiagramIndex>
MEqq2gZ2ffNason::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if ( diags[i]->id() == -1 ) sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -2 ) sel.insert(meInfo()[1], i);
  }
  return sel;
}

Selector<const ColourLines *>
MEqq2gZ2ffNason::colourGeometries(tcDiagPtr) const {
  static const ColourLines c1("1 -2");
  static const ColourLines c2("1 -2,4 -5");
  Selector<const ColourLines *> sel;
  if(abs(mePartonData()[2]->id())<=6) sel.insert(1.0, &c2);
  else                                sel.insert(1.0, &c1);
  return sel;
}

void MEqq2gZ2ffNason::persistentOutput(PersistentOStream & os) const {
  os << _maxflavour << _gammaZ << _process << _theFFZVertex << _theFFPVertex 
     << _gamma << _z0 << _contrib << _a << _p;
}

void MEqq2gZ2ffNason::persistentInput(PersistentIStream & is, int) { 
  is >> _maxflavour >> _gammaZ >> _process
     >> _theFFZVertex >> _theFFPVertex 
     >> _gamma >> _z0 >> _contrib >> _a >> _p; 
}

ClassDescription<MEqq2gZ2ffNason> MEqq2gZ2ffNason::initMEqq2gZ2ffNason;
// Definition of the static class description member.

void MEqq2gZ2ffNason::Init() {

  static ClassDocumentation<MEqq2gZ2ffNason> documentation
    ("The MEqq2gZ2ffNason class implements the matrix element for"
     "q qbar to Standard Model fermions via Z and photon exchange using"
     " helicity amplitude techniques");

  static Parameter<MEqq2gZ2ffNason,unsigned int> interfaceMaxFlavour
    ("MaxFlavour",
     "The heaviest incoming quark flavour this matrix element is allowed to handle",
     &MEqq2gZ2ffNason::_maxflavour, 5, 1, 6,
     false, false, Interface::limited);

  static Switch<MEqq2gZ2ffNason,unsigned int> interfaceGammaZ
    ("GammaZ",
     "Which terms to include",
     &MEqq2gZ2ffNason::_gammaZ, 0, false, false);
  static SwitchOption interfaceGammaZAll
    (interfaceGammaZ,
     "All",
     "Include both gamma and Z terms",
     0);
  static SwitchOption interfaceGammaZGamma
    (interfaceGammaZ,
     "Gamma",
     "Only include the photon",
     1);
  static SwitchOption interfaceGammaZZ
    (interfaceGammaZ,
     "Z",
     "Only include the Z",
     2);

  static Switch<MEqq2gZ2ffNason,unsigned int> interfaceProcess
    ("Process",
     "Which process to included",
     &MEqq2gZ2ffNason::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all SM fermions as outgoing particles",
     0);
  static SwitchOption interfaceProcessQuarks
    (interfaceProcess,
     "Quarks",
     "All include the quarks as outgoing particles",
     1);
  static SwitchOption interfaceProcessLeptons
    (interfaceProcess,
     "Leptons",
     "Only include the leptons as outgoing particles",
     2);
  static SwitchOption interfaceProcessChargedLeptons
    (interfaceProcess,
     "ChargedLeptons",
     "Only include the charged leptons as outgoing particles",
     3);
  static SwitchOption interfaceProcessNeutrinos
    (interfaceProcess,
     "Neutrinos",
     "Only include the neutrinos as outgoing particles",
     4);
  static SwitchOption interfaceProcessElectron
    (interfaceProcess,
     "Electron",
     "Only include e+e- as outgoing particles",
     5);
  static SwitchOption interfaceProcessMuon
    (interfaceProcess,
     "Muon",
     "Only include mu+mu- as outgoing particles",
     6);
  static SwitchOption interfaceProcessTau
    (interfaceProcess,
     "Tau",
     "Only include tau+tau- as outgoing particles",
     7);
  static SwitchOption interfaceProcessNu_e
    (interfaceProcess,
     "Nu_e",
     "Only include nu_e ne_ebar as outgoing particles",
     8);
  static SwitchOption interfaceProcessnu_mu
    (interfaceProcess,
     "nu_mu",
     "Only include nu_mu nu_mubar as outgoing particles",
     9);
  static SwitchOption interfaceProcessnu_tau
    (interfaceProcess,
     "nu_tau",
     "Only include nu_tau nu_taubar as outgoing particles",
     10);
  static SwitchOption interfaceProcessDown
    (interfaceProcess,
     "Down",
     "Only include d dbar as outgoing particles",
     11);
  static SwitchOption interfaceProcessUp
    (interfaceProcess,
     "Up",
     "Only include u ubar as outgoing particles",
     12);
  static SwitchOption interfaceProcessStrange
    (interfaceProcess,
     "Strange",
     "Only include s sbar as outgoing particles",
     13);
  static SwitchOption interfaceProcessCharm
    (interfaceProcess,
     "Charm",
     "Only include c cbar as outgoing particles",
     14);
  static SwitchOption interfaceProcessBottom
    (interfaceProcess,
     "Bottom",
     "Only include b bbar as outgoing particles",
     15);
  static SwitchOption interfaceProcessTop
    (interfaceProcess,
     "Top",
     "Only include t tbar as outgoing particles",
     16);

  static Switch<MEqq2gZ2ffNason,unsigned int> interfaceContribution
    ("Contribution",
     "Which contributions to the cross section to include",
     &MEqq2gZ2ffNason::_contrib, 1, false, false);
  static SwitchOption interfaceContributionLeadingOrder
    (interfaceContribution,
     "LeadingOrder",
     "Just generate the leading order cross section",
     0);
  static SwitchOption interfaceContributionPositiveNLO
    (interfaceContribution,
     "PositiveNLO",
     "Generate the positive contribution to the full NLO cross section",
     1);
  static SwitchOption interfaceContributionNegativeNLO
    (interfaceContribution,
     "NegativeNLO",
     "Generate the negative contribution to the full NLO cross section",
     2);

  static Parameter<MEqq2gZ2ffNason,double> interfaceCorrectionCoefficient
    ("CorrectionCoefficient",
     "The magnitude pf the correction term to reduce the negative contribution",
     &MEqq2gZ2ffNason::_a, 0.5, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<MEqq2gZ2ffNason,double> interfaceCorrectionPower
    ("CorrectionPower",
     "The power of the correction term to reduce the negative contribution",
     &MEqq2gZ2ffNason::_p, 0.7, 0.0, 1.0,
     false, false, Interface::limited);
}

double MEqq2gZ2ffNason::qqbarME(vector<SpinorWaveFunction>    & fin ,
			   vector<SpinorBarWaveFunction> & ain ,
			   vector<SpinorBarWaveFunction> & fout,
			   vector<SpinorWaveFunction>    & aout,
			   bool calc) const {
  // scale
  //Energy2 mb2(scale());
  Energy2 mb2(sqr(getParticleData(ParticleID::Z0)->mass()));
  // matrix element to be stored
  ProductionMatrixElement menew(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1Half,PDT::Spin1Half);
  // which intermediates to include
  bool gamma=_gammaZ==0||_gammaZ==1;
  bool Z0   =_gammaZ==0||_gammaZ==2;
  // declare the variables we need
  unsigned int ihel1,ihel2,ohel1,ohel2;
  VectorWaveFunction inter[2],test;
  double me[3]={0.,0.,0.};
  Complex diag1,diag2;
  // sum over helicities to get the matrix element
  for(ihel1=0;ihel1<2;++ihel1) {
    for(ihel2=0;ihel2<2;++ihel2) {
      // intermediate for Z
      if(Z0)    inter[0]=_theFFZVertex->evaluate(mb2,1,_z0   ,fin[ihel1],ain[ihel2]);
      // intermediate for photon
      if(gamma) inter[1]=_theFFPVertex->evaluate(mb2,1,_gamma,fin[ihel1],ain[ihel2]);
      for(ohel1=0;ohel1<2;++ohel1) {
	for(ohel2=0;ohel2<2;++ohel2) {
	  // first the Z exchange diagram
	  diag1 = Z0 ?
	    _theFFZVertex->evaluate(mb2,aout[ohel2],fout[ohel1],inter[0]) : 0.;
	  // first the photon exchange diagram
	  diag2 = gamma ?
	    _theFFPVertex->evaluate(mb2,aout[ohel2],fout[ohel1],inter[1]) : 0.;
	  // add up squares of individual terms
	  me[1] += real(diag1*conj(diag1));
	  me[2] += real(diag2*conj(diag2));
	  // the full thing including interference
	  diag1 +=diag2;
	  me[0] += real(diag1*conj(diag1));
	  if(calc) menew(ihel1,ihel2,ohel1,ohel2) = diag1;
	}
      }
    }
  }
  // spin and colour factor
  double colspin=1./12.;
  if(abs(fout[0].id())<=6) colspin*=3.;
  // results
  // factor 12 from 4 helicity and 3 colour
  for(int ix=0;ix<3;++ix){me[ix]*=colspin;}
  DVector save;
  save.push_back(me[1]);
  save.push_back(me[2]);
  meInfo(save);
  if(calc) _me.reset(menew);
  // the last factor is the difference betweeon our and MCFM's
  // alpha at the Z mass
  return me[0]*NLOweight()*sqr(0.00754677226/0.00776334);
}

void MEqq2gZ2ffNason::constructVertex(tSubProPtr sub) {
  SpinfoPtr spin[4];
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);hard.push_back(sub->outgoing()[1]);
  // order of particles
  unsigned int order[4]={0,1,2,3};
  if(hard[0]->id()<0){order[0]=1;order[1]=0;}
  if(hard[2]->id()<0){order[2]=3;order[3]=2;}
  vector<SpinorWaveFunction>    fin,aout;
  vector<SpinorBarWaveFunction> ain,fout;
  SpinorWaveFunction(   fin ,hard[order[0]],incoming,false,true);
  SpinorBarWaveFunction(ain ,hard[order[1]],incoming,false,true);
  SpinorBarWaveFunction(fout,hard[order[2]],outgoing,true ,true);
  SpinorWaveFunction(   aout,hard[order[3]],outgoing,true ,true);
  qqbarME(fin,ain,fout,aout,true);
  // get the spin info objects
  for(unsigned int ix=0;ix<4;++ix)
    {spin[ix]=dynamic_ptr_cast<SpinfoPtr>(hard[order[ix]]->spinInfo());}
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix){spin[ix]->setProductionVertex(hardvertex);}
}

int MEqq2gZ2ffNason::nDim() const {
  return 3;
}

bool MEqq2gZ2ffNason::generateKinematics(const double * r) {
  _x=*(r+1);
  _v=*(r+2);
  //_x=UseRandom::rnd();
  //_v=UseRandom::rnd();
  return ME2to2Base::generateKinematics(r);
}

double MEqq2gZ2ffNason::NLOweight() const {
  // if leading order this is one
  if(_contrib==0) return 1.;
  using Constants::pi;
  // first extract xbarp and xbarm
  double xbp = lastX1();
  double xbm = lastX2();
  tcPDPtr in[2]={mePartonData()[0],mePartonData()[1]};
  Ptr<BeamParticleData>::transient_const_pointer beam[2]=
    {dynamic_ptr_cast<Ptr<BeamParticleData>::transient_const_pointer>(lastParticles().first->dataPtr()),
     dynamic_ptr_cast<Ptr<BeamParticleData>::transient_const_pointer>(lastParticles().second->dataPtr())};
  tcPDPtr gluon = getParticleData(ParticleID::g);
  if(mePartonData()[0]->id()<0) {
    swap(xbp    ,xbm    );
    swap(beam[0],beam[1]);
    swap(in[0]  ,in[1]  );
  }
  // alpha_S
  const double eps=1e-8;
  double aS = SM().alphaS(scale());
  double CF(4./3.),TR(0.5);
  // pieces with LO kinematics
  double lowgt = 1.+aS/pi*CF;
  // + radiation and collinear pieces
  double z  = xbp+(1.-xbp)*_x;
  double xp = xbp/z;
  double oldq = beam[0]->pdf()->xfx(beam[0],in[0],scale(),xbp)/xbp;
  if(oldq<eps) return 0.;
  double newq = beam[0]->pdf()->xfx(beam[0],in[0],scale(),xp )/xp;
  double newg = beam[0]->pdf()->xfx(beam[0],gluon,scale(),xp )/xp;
  double vt = (1.-z)*_v; 
  double zlog = log(sqr(1.-z)/z);
  double pwgt = 0.5*aS/pi*(1.-xbp)/z*
    (2.*CF*(vt-(1.-z))     *newq/oldq
     +  TR*(1.-z)*(vt+2.*z)*newg/oldq
     +CF*((2.*sqr(pi)/3.-5.)+(1.-z-(1.+z)*zlog)*newq/oldq+
	  (newq/oldq-1.)*2./(1.-z)*zlog)
     +TR*((sqr(z)+sqr(1.-z))*zlog+2.*z*(1.-z))*newg/oldq);
  if(isnan(pwgt)||isinf(pwgt)) cerr << "testing + weight nan\n";
  int xbin = int(_x*100.);
  int vbin = int(_v*100.);
  if(pwgt>0.) {
    _posxp[xbin] += pwgt;
    _posvp[vbin] += pwgt; 
  }
  else {
    _negxp[xbin] += (-pwgt);
    _negvp[vbin] += (-pwgt); 
  }
  // - radiation and collinear pieces
  z  = xbm+(1.-xbm)*_x;
  double xm = xbm/z;
  oldq = beam[1]->pdf()->xfx(beam[1],in[1],scale(),xbm)/xbm;
  if(oldq<eps) return 0.;
  newq = beam[1]->pdf()->xfx(beam[1],in[1],scale(),xm )/xm;
  newg = beam[1]->pdf()->xfx(beam[1],gluon,scale(),xm )/xm;
  vt = (1.-z)*_v;
  zlog = log(sqr(1.-z)/z);
  double nwgt = 0.5*aS/pi*(1.-xbm)/z*
    (2.*CF*(vt-(1.-z))     *newq/oldq
     +  TR*(1.-z)*(vt+2.*z)*newg/oldq
     +CF*((2.*sqr(pi)/3.-5.)+(1.-z-(1.+z)*zlog)*newq/oldq+
	  (newq/oldq-1.)*2./(1.-z)*zlog)
     +TR*((sqr(z)+sqr(1.-z))*zlog+2.*z*(1.-z))*newg/oldq);
  if(isnan(nwgt)||isinf(nwgt)) cerr << "testing - weight nan\n";
  if(nwgt>0.) {
    _posxn[xbin] += nwgt;
    _posvn[vbin] += nwgt; 
  }
  else {
    _negxn[xbin] += (-nwgt);
    _negvn[vbin] += (-nwgt); 
  }
  double wgt = lowgt+pwgt+nwgt;
  // trick to try and make less negative events
  if(_x>=eps) {
    wgt += _a*(1./pow(_x,_p)-(1.-pow(eps,1.-_p))/(1.-_p)/(1.-eps));
  }
  if(xbin>99||vbin>99) cerr << "testing hist error\n";
  if(wgt>0.) {
    _posx[xbin] += wgt;
    _posv[vbin] += wgt; 
  }
  else {
    _negx[xbin] += (-wgt);
    _negv[vbin] += (-wgt); 
  }
  if(isnan(wgt)||isinf(wgt)) {
    cerr << "testing infinite weightA " << lowgt+pwgt+nwgt << "\n";
    cerr << "testing infinite weightB " << wgt << "\n";
    cerr << "testing variables " << _x << " " << _v << "\n";
    exit(0);
  }
  return _contrib==1 ? max(0.,wgt) : max(0.,-wgt);
}

void MEqq2gZ2ffNason::dofinish() {
  ME2to2Base::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"weight for x distribution +ve weight events\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"x\"\n";
  double step=0.01;
  double x=-0.005;
  for(unsigned int ix=0;ix<_posx.size();++ix) {
    x+=step;
    outfile << x << "\t" << _posx[ix].mean() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"weight for x distribution -ve weight events\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"x\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_negx.size();++ix) {
    x+=step;
    outfile << x << "\t" << _negx[ix].mean() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"weight for v distribution +ve weight events\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"v\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_posx.size();++ix) {
    x+=step;
    outfile << x << "\t" << _posv[ix].mean() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"weight for v distribution -ve weight events\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"v\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_negx.size();++ix) {
    x+=step;
    outfile << x << "\t" << _negv[ix].mean() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"numbers of enteries for x distribution +ve weight events\"\n";
  outfile << "TITLE LEFT \"numbers of enteries\"\n";
  outfile << "TITLE BOTTOM \"x\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_posx.size();++ix) {
    x+=step;
    outfile << x << "\t" << _posx[ix].numberOfPoints() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"numbers of enteries for x distribution -ve weight events\"\n";
  outfile << "TITLE LEFT \"numbers of enteries\"\n";
  outfile << "TITLE BOTTOM \"x\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_negx.size();++ix) {
    x+=step;
    outfile << x << "\t" << _negx[ix].numberOfPoints() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"numbers of enteries for v distribution +ve weight events\"\n";
  outfile << "TITLE LEFT \"numbers of enteries\"\n";
  outfile << "TITLE BOTTOM \"v\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_posx.size();++ix) {
    x+=step;
    outfile << x << "\t" << _posv[ix].numberOfPoints() << "\n";
  }
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"weight for x distribution +ve weight events Pos emission\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"x\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_posxp.size();++ix) {
    x+=step;
    outfile << x << "\t" << _posxp[ix].mean() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"weight for x distribution -ve weight events Pos emission\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"x\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_negxp.size();++ix) {
    x+=step;
    outfile << x << "\t" << _negxp[ix].mean() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"weight for v distribution +ve weight events Pos emission\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"v\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_posxp.size();++ix) {
    x+=step;
    outfile << x << "\t" << _posvp[ix].mean() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"weight for v distribution -ve weight events Pos emission\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"v\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_negxp.size();++ix) {
    x+=step;
    outfile << x << "\t" << _negvp[ix].mean() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"numbers of enteries for x distribution +ve weight events Pos emission\"\n";
  outfile << "TITLE LEFT \"numbers of enteries\"\n";
  outfile << "TITLE BOTTOM \"x\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_posxp.size();++ix) {
    x+=step;
    outfile << x << "\t" << _posxp[ix].numberOfPoints() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"numbers of enteries for x distribution -ve weight events Pos emission\"\n";
  outfile << "TITLE LEFT \"numbers of enteries\"\n";
  outfile << "TITLE BOTTOM \"x\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_negxp.size();++ix) {
    x+=step;
    outfile << x << "\t" << _negxp[ix].numberOfPoints() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"numbers of enteries for v distribution +ve weight events Pos emission\"\n";
  outfile << "TITLE LEFT \"numbers of enteries\"\n";
  outfile << "TITLE BOTTOM \"v\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_posxp.size();++ix) {
    x+=step;
    outfile << x << "\t" << _posvp[ix].numberOfPoints() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"numbers of enteries for v distribution -ve weight events Pos emission\"\n";
  outfile << "TITLE LEFT \"numbers of enteries\"\n";
  outfile << "TITLE BOTTOM \"v\"\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"weight for x distribution +ve weight events Neg emission\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"x\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_posxn.size();++ix) {
    x+=step;
    outfile << x << "\t" << _posxn[ix].mean() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"weight for x distribution -ve weight events Neg emission\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"x\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_negxn.size();++ix) {
    x+=step;
    outfile << x << "\t" << _negxn[ix].mean() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"weight for v distribution +ve weight events Neg emission\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"v\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_posxn.size();++ix) {
    x+=step;
    outfile << x << "\t" << _posvn[ix].mean() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"weight for v distribution -ve weight events Neg emission\"\n";
  outfile << "TITLE LEFT \"weight\"\n";
  outfile << "TITLE BOTTOM \"v\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_negxn.size();++ix) {
    x+=step;
    outfile << x << "\t" << _negvn[ix].mean() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"numbers of enteries for x distribution +ve weight events Neg emission\"\n";
  outfile << "TITLE LEFT \"numbers of enteries\"\n";
  outfile << "TITLE BOTTOM \"x\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_posxn.size();++ix) {
    x+=step;
    outfile << x << "\t" << _posxn[ix].numberOfPoints() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"numbers of enteries for x distribution -ve weight events Neg emission\"\n";
  outfile << "TITLE LEFT \"numbers of enteries\"\n";
  outfile << "TITLE BOTTOM \"x\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_negxn.size();++ix) {
    x+=step;
    outfile << x << "\t" << _negxn[ix].numberOfPoints() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"numbers of enteries for v distribution +ve weight events Neg emission\"\n";
  outfile << "TITLE LEFT \"numbers of enteries\"\n";
  outfile << "TITLE BOTTOM \"v\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_posxn.size();++ix) {
    x+=step;
    outfile << x << "\t" << _posvn[ix].numberOfPoints() << "\n";
  }
  outfile << "HIST\n";
  outfile << "NEW FRAME\n";
  outfile << "TITLE TOP \"numbers of enteries for v distribution -ve weight events Neg emission\"\n";
  outfile << "TITLE LEFT \"numbers of enteries\"\n";
  outfile << "TITLE BOTTOM \"v\"\n";
  x=-0.005;
  for(unsigned int ix=0;ix<_negxn.size();++ix) {
    x+=step;
    outfile << x << "\t" << _negvn[ix].numberOfPoints() << "\n";
  }
  outfile << "HIST\n";
  outfile.close();
}
