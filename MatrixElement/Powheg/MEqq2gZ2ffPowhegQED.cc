// -*- C++ -*-
  
// This is the implementation of the non-inlined, non-templated member
// functions of the MEqq2gZ2ffPowhegQED class.
//

#include "MEqq2gZ2ffPowhegQED.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "Herwig++/Utilities/Maths.h"
#include "Herwig++/MatrixElement/HardVertex.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include <numeric>


using namespace Herwig;

/**
 *  Typedef for BeamParticleData
 */
typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr;

MEqq2gZ2ffPowhegQED::MEqq2gZ2ffPowhegQED() 
  : corrections_(1), QEDContributions_(0), incomingPhotons_(true),
    DipoleSum_(0),
    contrib_(1), power_(0.6), zPow_(0.5), yPow_(0.9),
    alphaS_(0.118), alphaEM_(1./137.), fixedCouplings_(false),
    supressionFunction_(0), 
    lambda2_(10000.*GeV2),QCDRadiationType_(-1),
    preqqbarqQCD_(10.), preqqbarqbarQCD_(10.),
    preqgQCD_(10.),pregqbarQCD_(10.),
    preqqbarqQED_(50.), preqqbarqbarQED_(50.),
    preqgQED_(10.),pregqbarQED_(10.), preFFQED_(6.),
    preIFQED_(10.),
    minpTQCD_(2.*GeV),
    minpTQEDII_(1.*GeV),minpTQEDIF_(1.*GeV),minpTQEDFF_(1.*GeV),
    process_(2), maxFlavour_(5) {
  vector<unsigned int> mopt(2,1);
  massOption(mopt);
}

unsigned int MEqq2gZ2ffPowhegQED::orderInAlphaS() const {
  return 0;
}

unsigned int MEqq2gZ2ffPowhegQED::orderInAlphaEW() const {
  return 2;
}

IBPtr MEqq2gZ2ffPowhegQED::clone() const {
  return new_ptr(*this);
}

IBPtr MEqq2gZ2ffPowhegQED::fullclone() const {
  return new_ptr(*this);
}

Energy2 MEqq2gZ2ffPowhegQED::scale() const {
  return sqr(Z0_->mass());
  //return sHat();
}

int MEqq2gZ2ffPowhegQED::nDim() const {
  return HwMEBase::nDim() + ( contrib_>=1 && contrib_<=3 ? 3 : 0 );
}

bool MEqq2gZ2ffPowhegQED::generateKinematics(const double * r) {
  if(contrib_>=1&&contrib_<=3) {
    zTilde_ = r[nDim()-1];
    vTilde_ = r[nDim()-2];
    phi_    = Constants::twopi*r[nDim()-3];
  }
  jacobian(1.0);
  return HwMEBase::generateKinematics(r);
}

CrossSection MEqq2gZ2ffPowhegQED::dSigHatDR() const {
  // old technique
  CrossSection preFactor = 
    jacobian()/(16.0*sqr(Constants::pi)*sHat())*sqr(hbarc);
  loME_ = me2();
  if(contrib_<4) return NLOWeight()*preFactor;
  // folding technique to ensure positive
  double wgt(0.);
  unsigned int ntry(0);
  do {
    // radiative variables
    zTilde_ = UseRandom::rnd();
    vTilde_ = UseRandom::rnd();
    phi_    = Constants::twopi*UseRandom::rnd();
    wgt += NLOWeight();
    ++ntry;
  }
  while (wgt<0.&&ntry<100);
  if(wgt<0.) return ZERO;
  return wgt*preFactor/double(ntry);
}

double MEqq2gZ2ffPowhegQED::loME(const cPDVector & particles,
				 const vector<Lorentz5Momentum> & momenta,
				 bool first) const {
  vector<SpinorWaveFunction>    fin,aout;
  vector<SpinorBarWaveFunction> ain,fout;
  SpinorWaveFunction       q(momenta[0],particles[0],incoming);
  SpinorBarWaveFunction qbar(momenta[1],particles[1],incoming);
  SpinorBarWaveFunction    f(momenta[2],particles[2],outgoing);
  SpinorWaveFunction    fbar(momenta[3],particles[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    q.reset(ix)   ; fin.push_back(q);
    qbar.reset(ix); ain.push_back(qbar);
    f.reset(ix)   ;fout.push_back(f);
    fbar.reset(ix);aout.push_back(fbar);
  }
  return qqbarME(fin,ain,fout,aout,first);
}

double MEqq2gZ2ffPowhegQED::qqbarME(vector<SpinorWaveFunction>    & fin ,
				    vector<SpinorBarWaveFunction> & ain ,
				    vector<SpinorBarWaveFunction> & fout,
				    vector<SpinorWaveFunction>    & aout,
				    bool first) const {
  // scale for the process
  const Energy2 q2(scale());
  ProductionMatrixElement menew(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1Half,PDT::Spin1Half);
  VectorWaveFunction inter[2];
  double me[3]={0.,0.,0.};
  // sum over helicities to get the matrix element
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      // intermediate for Z
      inter[0]=oneLoopVertex_->evaluate(q2,1,Z0_   ,fin[ihel1],ain[ihel2]);
      // intermediate for photon
      inter[1]=oneLoopVertex_->evaluate(q2,1,gamma_,fin[ihel1],ain[ihel2]);
      for(unsigned int ohel1=0;ohel1<2;++ohel1) {
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
 	  // first the Z exchange diagram
 	  Complex diag1 = oneLoopVertex_->
	    evaluate(q2,aout[ohel2],fout[ohel1],inter[0]);
	  // first the photon exchange diagram
	  Complex diag2 =  oneLoopVertex_->
	    evaluate(q2,aout[ohel2],fout[ohel1],inter[1]);
	  // add up squares of individual terms
	  me[1] += norm(diag1);
	  me[2] += norm(diag2);
	  // the full thing including interference
 	  diag1 += diag2;
	  me[0] += norm(diag1);
	  menew(ihel1,ihel2,ohel1,ohel2) = diag1;
	}
      }
    }
  }
  // spin and colour factor
  double colspin=1./12.;
  if(abs(fout[0].id())<=6) colspin*=3.;
  for(int ix=0;ix<3;++ix) 
    me[ix]*=colspin;
  if(first) {
    DVector save;
    save.push_back(me[1]);
    save.push_back(me[2]);
    meInfo(save);
    me_.reset(menew);
  }
  return me[0];
}

Selector<const ColourLines *>
MEqq2gZ2ffPowhegQED::colourGeometries(tcDiagPtr) const {
  static const ColourLines c1("1 -2");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &c1);
  return sel;
}

void MEqq2gZ2ffPowhegQED::getDiagrams() const {
  // loop over the processes we need
  for ( int ix=11; ix<17; ++ix ) {
    // is it a valid lepton process
    bool lepton= 
      ( process_==2              || (process_==3 && ix%2==1) || 
	(process_==4 && ix%2==0) || (ix%2==0 && (ix-10)/2==process_-7) ||
	(ix%2==1 && (ix-9)/2 ==process_-4));
    // if not a valid process continue
    if(!lepton) continue;
    tcPDPtr lm = getParticleData(ix);
    tcPDPtr lp = lm->CC();
    for(int i = 1; i <= maxFlavour_; ++i) {
      tcPDPtr q  = getParticleData(i);
      tcPDPtr qb = q->CC();
      add(new_ptr((Tree2toNDiagram(2), q, qb, 1, Z0_   , 3, lm, 3, lp, -1)));
      add(new_ptr((Tree2toNDiagram(2), q, qb, 1, gamma_, 3, lm, 3, lp, -2)));
    }
  }
}

Selector<MEBase::DiagramIndex>
MEqq2gZ2ffPowhegQED::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if ( diags[i]->id() == -1 ) sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -2 ) sel.insert(meInfo()[1], i);
  }
  return sel;
}

double MEqq2gZ2ffPowhegQED::me2() const {
  // clear looptools cache if needed
  if(contrib_!=0&&(corrections_==3||corrections_==5))
    oneLoopVertex_->clearCache();
  // compute the spinors
  vector<SpinorWaveFunction>     fin,aout;
  vector<SpinorBarWaveFunction>  ain,fout;
  SpinorWaveFunction       q(meMomenta()[0],mePartonData()[0],incoming);
  SpinorBarWaveFunction qbar(meMomenta()[1],mePartonData()[1],incoming);
  SpinorBarWaveFunction    f(meMomenta()[2],mePartonData()[2],outgoing);
  SpinorWaveFunction    fbar(meMomenta()[3],mePartonData()[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    q.reset(ix)   ;
    fin .push_back(q);
    qbar.reset(ix);
    ain .push_back(qbar);
    f   .reset(ix);
    fout.push_back(f);
    fbar.reset(ix);
    aout.push_back(fbar);
  }
  // scale for the process
  const Energy2 q2(scale());
  // wavefunctions and box corrections
  // for EW correction if needed
  VectorWaveFunction FS_Z[2][2],FS_A[2][2];
  VectorWaveFunction IS_Z[2][2],IS_A[2][2];
  vector<vector<complex<InvEnergy2 > > > boxCoeff;
  if(contrib_!=0&&(corrections_==3||corrections_==5)) {
    // compute final-state vertex correction
    oneLoopVertex_->setOrder(1);
    for(unsigned int ohel1=0;ohel1<2;++ohel1) {
      for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	FS_Z[ohel1][ohel2] = oneLoopVertex_->
	  evaluate(q2,1,Z0_   ,aout[ohel2],fout[ohel1]);
	FS_A[ohel1][ohel2] = oneLoopVertex_->
	  evaluate(q2,1,gamma_,aout[ohel2],fout[ohel1]);
      }
    }
    // compute initial-state vertex correction
    for(unsigned int ihel1=0;ihel1<2;++ihel1) {
      for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	IS_Z[ihel1][ihel2] = oneLoopVertex_->
	  evaluate(q2,1,Z0_   ,fin[ihel1],ain[ihel2]);
	IS_A[ihel1][ihel2] = oneLoopVertex_->
	  evaluate(q2,1,gamma_,fin[ihel1],ain[ihel2]);
      }
    }
    oneLoopVertex_->setOrder(0);
    // coefficients for the box diagrams
    boxCoeff = oneLoopVertex_->
      neutralCurrentFBox(mePartonData()[0],mePartonData()[1],
			 mePartonData()[2],mePartonData()[3],
			 sHat(),tHat(),uHat());
  }
  // momentum difference for genuine NLO QED FS structure
  LorentzPolarizationVector momDiff = 
    (rescaledMomenta()[2]-rescaledMomenta()[3])/2./
    (rescaledMomenta()[2].mass()+rescaledMomenta()[3].mass());
  // storage of ME
  ProductionMatrixElement menew(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1Half,PDT::Spin1Half);
  // intermediate wavefunctions
  VectorWaveFunction inter_LO_AA ,inter_LO_ZZ ;
  VectorWaveFunction inter_NLO_AA,inter_NLO_AZ,
    inter_NLO_ZA,inter_NLO_ZZ;
  // sum over helicities to get the matrix element
  double loSum(0.),EWSum(0.),QEDSum(0.);
  double  dsum[2]={0.,0.};
  // incoming helicities
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      // LO intermediates
      // intermediate for Z
      inter_LO_ZZ  = oneLoopVertex_->evaluate(q2,1,Z0_   ,fin[ihel1],ain[ihel2]);
      // intermediate for photon
      inter_LO_AA  = oneLoopVertex_->evaluate(q2,1,gamma_,fin[ihel1],ain[ihel2]);
      // stuff for EW correction if needed
      LorentzPolarizationVectorE leftQuark,rightQuark;
      if(contrib_!=0&&(corrections_==3||corrections_==5)) {
	// intermediates including self energy corrections 
	inter_NLO_AA = oneLoopVertex_->selfEnergyCorrection(gamma_,inter_LO_AA);
	inter_NLO_AZ = oneLoopVertex_->selfEnergyCorrection(Z0_   ,inter_LO_AA);
	inter_NLO_ZA = oneLoopVertex_->selfEnergyCorrection(gamma_,inter_LO_ZZ);
	inter_NLO_ZZ = oneLoopVertex_->selfEnergyCorrection(Z0_   ,inter_LO_ZZ);
	// currents for the box diagrams
	leftQuark  = fin[ihel1].dimensionedWave().
	   leftCurrent(ain[ihel2].dimensionedWave());
	rightQuark = fin[ihel1].dimensionedWave().
	  rightCurrent(ain[ihel2].dimensionedWave());
      }
      // scalars for QED correction if needed
      Complex scalar1(0.),scalar2(0.);
      if(contrib_!=0&&(corrections_==2||corrections_==5)) {
	scalar1 = inter_LO_ZZ.wave().dot(momDiff);
	scalar2 = inter_LO_AA.wave().dot(momDiff);
      }
      // outgoing helicities
      for(unsigned int ohel1=0;ohel1<2;++ohel1) {
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  // LO piece
 	  // first the Z exchange diagram
 	  Complex diag1 = oneLoopVertex_->
	    evaluate(q2,aout[ohel2],fout[ohel1],inter_LO_ZZ);
	  // first the photon exchange diagram
	  Complex diag2 =  oneLoopVertex_->
	    evaluate(q2,aout[ohel2],fout[ohel1],inter_LO_AA);
	  Complex lo = diag1 + diag2;
	  loSum += norm(lo);
 	  // add up squares of individual terms
 	  dsum[0] += norm(diag1);
 	  dsum[1] += norm(diag2);
	  menew(ihel1,ihel2,ohel1,ohel2) = lo;
	  // genuine EW piece
	  if(contrib_!=0&&(corrections_==3||corrections_==5)) {
	    // self energy piece
	    Complex diag3  = oneLoopVertex_->
	      evaluate(q2,aout[ohel2],fout[ohel1],inter_NLO_AA);
	    Complex diag4  = oneLoopVertex_->
	      evaluate(q2,aout[ohel2],fout[ohel1],inter_NLO_AZ);
	    Complex diag5  = oneLoopVertex_->
	      evaluate(q2,aout[ohel2],fout[ohel1],inter_NLO_ZA);
	    Complex diag6  = oneLoopVertex_->
	      evaluate(q2,aout[ohel2],fout[ohel1],inter_NLO_ZZ);
	    Complex self   = diag3 + diag4 + diag5 + diag6;
	    // vertex corrections
	    Complex diag7  = oneLoopVertex_->
	      evaluate(q2,aout[ohel2],fout[ohel1],IS_Z[ihel1][ihel2]);
	    Complex diag8  = oneLoopVertex_->
	      evaluate(q2,aout[ohel2],fout[ohel1],IS_A[ihel1][ihel2]);
	    Complex diag9  = oneLoopVertex_->
	      evaluate(q2,fin [ihel1],ain [ihel2],FS_Z[ohel1][ohel2]);
	    Complex diag10 = oneLoopVertex_->
	      evaluate(q2,fin [ihel1],ain [ihel2],FS_A[ohel1][ohel2]);
	    Complex vertex = diag7 + diag8 + diag9 + diag10;
	    // box diagrams
	    // currents
	    LorentzPolarizationVectorE leftLepton  = 
	      aout[ohel2].dimensionedWave().leftCurrent (fout[ohel1].dimensionedWave());
	    LorentzPolarizationVectorE rightLepton = 
	      aout[ohel2].dimensionedWave().rightCurrent(fout[ohel1].dimensionedWave());
	    Complex box = -Complex(0.,1.)*
	      ( leftQuark.dot( leftLepton)*boxCoeff[0][0] + 
		leftQuark.dot(rightLepton)*boxCoeff[0][1] +
		rightQuark.dot( leftLepton)*boxCoeff[1][0] + 
		rightQuark.dot(rightLepton)*boxCoeff[1][1]);
	    // sum of the corrections
	    Complex loop = vertex + self + box; 
	    EWSum += real(loop*conj(lo)+lo*conj(loop));
	  }
	  // QED FS corrections
	  if(contrib_!=0&&(corrections_==2||corrections_==5)) {
	    LorentzPolarizationVector left  = 
	      aout[ohel2].wave(). leftCurrent(fout[ohel1].wave());
	    LorentzPolarizationVector right = 
	      aout[ohel2].wave().rightCurrent(fout[ohel1].wave());
	    Complex scalar = 
	      aout[ohel2].wave().scalar(fout[ohel1].wave());
	    // nlo specific pieces
	    Complex diag3 =
	      Complex(0.,1.)*oneLoopVertex_->norm()*
	      (oneLoopVertex_->right()*( left.dot(inter_LO_ZZ.wave())) +
	       oneLoopVertex_-> left()*(right.dot(inter_LO_ZZ.wave())) -
	       ( oneLoopVertex_-> left()+oneLoopVertex_->right())*scalar1*scalar);
	    diag3 += Complex(0.,1.)*oneLoopVertex_->norm()*
	      (oneLoopVertex_->right()*( left.dot(inter_LO_AA.wave())) +
	       oneLoopVertex_-> left()*(right.dot(inter_LO_AA.wave())) -
	       ( oneLoopVertex_-> left()+oneLoopVertex_->right())*scalar2*scalar);
	    // nlo piece
	    QEDSum += real(diag1*conj(diag3) + diag3*conj(diag1));
	  }
	}
      }
    }
  }
  // spin and colour factor
  double colSpin=1./12.;
  if(abs(mePartonData()[2]->id())<=6) colSpin *= 3.;
  loSum  *= colSpin;
  EWSum  *= colSpin;
  QEDSum *= colSpin;
  for(unsigned int ix=0;ix<2;++ix) dsum[ix] *= colSpin;
  // test of EW correction (if needed)
//   cerr << "testing in EW " << loME_ << " " << loSum << " " << EWSum << "\n";
//   oneLoopVertex_->neutralCurrentEWME(mePartonData()[0],mePartonData()[1],
// 				   mePartonData()[2],mePartonData()[3],
// 				   sHat(),tHat(),uHat());
  // save the stuff for diagram selection
  DVector save;
  save.push_back(dsum[0]);
  save.push_back(dsum[1]);
  meInfo(save);
  me_.reset(menew);
  // save the higher order pieces
  f2term_ = QEDSum;
  EWterm_ = EWSum/loSum;
  // return LO result
  return loSum;
}

void MEqq2gZ2ffPowhegQED::persistentOutput(PersistentOStream & os) const {
  os << corrections_ << incomingPhotons_ << contrib_ << power_ << gluon_ 
     << fixedCouplings_ << alphaS_ << alphaEM_ << yPow_ << zPow_
     << supressionFunction_ << ounit(lambda2_,GeV2)
     << preqqbarqQCD_ << preqqbarqbarQCD_ << preqgQCD_ << pregqbarQCD_
     << preqqbarqQED_ << preqqbarqbarQED_ << preqgQED_ << pregqbarQED_
     << preFFQED_ << preIFQED_ << prefactorQCD_ << prefactorQED_ 
     << ounit(minpTQCD_,GeV) 
     << ounit(minpTQEDII_,GeV) << ounit(minpTQEDIF_,GeV) << ounit(minpTQEDFF_,GeV) 
     << alphaQCD_ << alphaQED_
     << FFGVertex_ << oneLoopVertex_ << QCDRadiationType_
     << Z0_ << gamma_ << process_ << maxFlavour_ << QEDContributions_ << DipoleSum_;
}

void MEqq2gZ2ffPowhegQED::persistentInput(PersistentIStream & is, int) {
  is >> corrections_ >> incomingPhotons_ >> contrib_ >> power_ >> gluon_ 
     >> fixedCouplings_ >> alphaS_ >> alphaEM_ >> yPow_ >> zPow_
     >> supressionFunction_ >> iunit(lambda2_,GeV2)
     >> preqqbarqQCD_ >> preqqbarqbarQCD_ >> preqgQCD_ >> pregqbarQCD_
     >> preqqbarqQED_ >> preqqbarqbarQED_ >> preqgQED_ >> pregqbarQED_
     >> preFFQED_ >> preIFQED_ >> prefactorQCD_ >> prefactorQED_ 
     >> iunit(minpTQCD_,GeV)
     >> iunit(minpTQEDII_,GeV) >> iunit(minpTQEDIF_,GeV) >> iunit(minpTQEDFF_,GeV) 
     >> alphaQCD_ >> alphaQED_
     >> FFGVertex_ >> oneLoopVertex_ >> QCDRadiationType_
     >> Z0_ >> gamma_ >> process_ >> maxFlavour_ >> QEDContributions_ >> DipoleSum_;
}

void MEqq2gZ2ffPowhegQED::doinit() {
  HwMEBase::doinit();
  // get the ParticleData objects
  Z0_    = getParticleData(ThePEG::ParticleID::Z0);
  gamma_ = getParticleData(ThePEG::ParticleID::gamma);
  gluon_ = getParticleData(ParticleID::g);
  // prefactors for overestimate for real emission
  // QCD
  prefactorQCD_.push_back(preqqbarqQCD_);
  prefactorQCD_.push_back(preqqbarqbarQCD_);
  prefactorQCD_.push_back(preqgQCD_);
  prefactorQCD_.push_back(pregqbarQCD_);
  // QED
  prefactorQED_.push_back(preqqbarqQED_);
  prefactorQED_.push_back(preqqbarqbarQED_);
  prefactorQED_.push_back(preqgQED_);
  prefactorQED_.push_back(pregqbarQED_);
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm)
    throw InitException() << "Must be the Herwig++ SusyBase class in "
			  << "MEqq2gZ2ffPowhegQED::doinit" 
			  << Exception::abortnow;
  // extract the vertices
  FFGVertex_ = hwsm->vertexFFG();
}

void MEqq2gZ2ffPowhegQED::doinitrun() {
  oneLoopVertex_->initrun();
  HwMEBase::doinitrun();
}

ClassDescription<MEqq2gZ2ffPowhegQED> MEqq2gZ2ffPowhegQED::initMEqq2gZ2ffPowhegQED;
// Definition of the static class description member.

void MEqq2gZ2ffPowhegQED::Init() {

  static ClassDocumentation<MEqq2gZ2ffPowhegQED> documentation
    ("The MEqq2gZ2ffPowhegQED class implements the strong and QED corrections"
     " to the Drell-Yan process");

  static Switch<MEqq2gZ2ffPowhegQED,unsigned int> interfaceCorrections
    ("Corrections",
     "Which corrections to include",
     &MEqq2gZ2ffPowhegQED::corrections_, 1, false, false);
  static SwitchOption interfaceCorrectionsQCD
    (interfaceCorrections,
     "QCD",
     "Only include the QCD corrections",
     1);
  static SwitchOption interfaceCorrectionsQED
    (interfaceCorrections,
     "QED",
     "Only include the QED corrections",
     2);
  static SwitchOption interfaceCorrectionsEW
    (interfaceCorrections,
     "EW",
     "Only include the EW corrections",
     3);
  static SwitchOption interfaceCorrectionsQCDandQED
    (interfaceCorrections,
     "QCDandQED",
     "Include both QED and QCD corrections",
     4);
  static SwitchOption interfaceCorrectionsAll
    (interfaceCorrections,
     "All",
     "Include QED, QCD and EW corrections",
     5);

  static Switch<MEqq2gZ2ffPowhegQED,bool> interfaceIncomingPhotons
    ("IncomingPhotons",
     "Whether or not to include incoming photons",
     &MEqq2gZ2ffPowhegQED::incomingPhotons_, true, false, false);
  static SwitchOption interfaceIncomingPhotonsYes
    (interfaceIncomingPhotons,
     "Yes",
     "Include them",
     true);
  static SwitchOption interfaceIncomingPhotonsNo
    (interfaceIncomingPhotons,
     "No",
     "Don't include them",
     false);

  static Switch<MEqq2gZ2ffPowhegQED,unsigned int> interfaceContribution
    ("Contribution",
     "Which contributions to the cross section to include",
     &MEqq2gZ2ffPowhegQED::contrib_, 1, false, false);
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
  static SwitchOption interfaceContributionReal
    (interfaceContribution,
     "Real",
     "Generate the high pT real emission only",
     3);
  static SwitchOption interfaceContributionTotalNLO
    (interfaceContribution,
     "TotalNLO",
     "Generate the full NLO cross section using folding",
     4);

  static Parameter<MEqq2gZ2ffPowhegQED,double> interfaceSamplingPower
    ("SamplingPower",
     "Power for the sampling of xp",
     &MEqq2gZ2ffPowhegQED::power_, 0.6, 0.0, 1.,
     false, false, Interface::limited);

  static Switch<MEqq2gZ2ffPowhegQED,bool> interfaceFixedAlphaS
    ("FixedAlpha",
     "Use a fixed value of alpha_S and alpha_EM",
     &MEqq2gZ2ffPowhegQED::fixedCouplings_, false, false, false);
  static SwitchOption interfaceFixedAlphaSYes
    (interfaceFixedAlphaS,
     "Yes",
     "Use fixed alpha_S and alpha_EM",
     true);
  static SwitchOption interfaceFixedAlphaSNo
    (interfaceFixedAlphaS,
     "No",
     "Use running alpha_S and alpha_EM",
     false);

  static Parameter<MEqq2gZ2ffPowhegQED,double> interfaceAlphaS
    ("AlphaS",
     "The fixed value of alpha_S to use",
     &MEqq2gZ2ffPowhegQED::alphaS_, 0., 0., 1.,
     false, false, Interface::limited);

  static Parameter<MEqq2gZ2ffPowhegQED,double> interfaceAlphaEM
    ("AlphaEM",
     "The fixed value of alpha_EM to use",
     &MEqq2gZ2ffPowhegQED::alphaEM_, 1./137., 0., 1.,
     false, false, Interface::limited);

  static Switch<MEqq2gZ2ffPowhegQED,unsigned int> interfaceSupressionFunction
    ("SupressionFunction",
     "Choice of the supression function",
     &MEqq2gZ2ffPowhegQED::supressionFunction_, 0, false, false);
  static SwitchOption interfaceSupressionFunctionNone
    (interfaceSupressionFunction,
     "None",
     "Default POWHEG approach",
     0);
  static SwitchOption interfaceSupressionFunctionThetaFunction
    (interfaceSupressionFunction,
     "ThetaFunction",
     "Use theta functions at scale Lambda",
     1);
  static SwitchOption interfaceSupressionFunctionSmooth
    (interfaceSupressionFunction,
     "Smooth",
     "Supress high pT by pt^2/(pt^2+lambda^2)",
     2);

  static Parameter<MEqq2gZ2ffPowhegQED,Energy2> interfaceSupressionScale
    ("SupressionScale",
     "The square of the scale for the supression function",
     &MEqq2gZ2ffPowhegQED::lambda2_, GeV2, 10000.0*GeV2, 0.0*GeV2, 0*GeV2,
     false, false, Interface::lowerlim);

  static Parameter<MEqq2gZ2ffPowhegQED,double> interfaceQQbarQPreFactorQCD
    ("QQbarQPreFactorQCD",
     "Prefactor for the sampling on qqbar -> X g with radiation from q",
     &MEqq2gZ2ffPowhegQED::preqqbarqQCD_, 20.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<MEqq2gZ2ffPowhegQED,double> interfaceQQbarQbarPreFactorQCD
    ("QQbarQbarPreFactorQCD",
     "Prefactor for the sampling on qqbar -> X g with radiation from qbar",
     &MEqq2gZ2ffPowhegQED::preqqbarqbarQCD_, 20.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<MEqq2gZ2ffPowhegQED,double> interfaceQGPreFactorQCD
    ("QGPreFactorQCD",
     "The prefactor for the qg->Xq channel",
     &MEqq2gZ2ffPowhegQED::preqgQCD_, 20.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<MEqq2gZ2ffPowhegQED,double> interfaceQbarGPreFactorQCD
    ("QbarGPreFactorQCD",
     "The prefactor for the qbarg->Xqbar channel",
     &MEqq2gZ2ffPowhegQED::pregqbarQCD_, 20.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<MEqq2gZ2ffPowhegQED,double> interfaceQQbarQPreFactorQED
    ("QQbarQPreFactorQED",
     "Prefactor for the sampling on qqbar -> X g with radiation from q",
     &MEqq2gZ2ffPowhegQED::preqqbarqQED_, 20.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<MEqq2gZ2ffPowhegQED,double> interfaceQQbarQbarPreFactorQED
    ("QQbarQbarPreFactorQED",
     "Prefactor for the sampling on qqbar -> X g with radiation from qbar",
     &MEqq2gZ2ffPowhegQED::preqqbarqbarQED_, 20.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<MEqq2gZ2ffPowhegQED,double> interfaceQGPreFactorQED
    ("QGPreFactorQED",
     "The prefactor for the qg->Xq channel",
     &MEqq2gZ2ffPowhegQED::preqgQED_, 20.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<MEqq2gZ2ffPowhegQED,double> interfaceQbarGPreFactorQED
    ("QbarGPreFactorQED",
     "The prefactor for the qbarg->Xqbar channel",
     &MEqq2gZ2ffPowhegQED::pregqbarQED_, 20.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<MEqq2gZ2ffPowhegQED,double> interfaceFinalStatePreFactorQED
    ("FinalStatePreFactorQED",
     "The prefactor for final-state QED radiation",
     &MEqq2gZ2ffPowhegQED::preFFQED_, 6.0, 0.0, 100.0,
     false, false, Interface::limited);

  static Parameter<MEqq2gZ2ffPowhegQED,double> interfaceInitialFinalStatePreFactorQED
    ("InitialFinalStatePreFactorQED",
     "The prefactor for inital- and final-state QED interference",
     &MEqq2gZ2ffPowhegQED::preIFQED_, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static Parameter<MEqq2gZ2ffPowhegQED,Energy> interfaceMinimumpTQCD
    ("MinimumpTQCD",
     "The minimum pT for the hard QCD emission",
     &MEqq2gZ2ffPowhegQED::minpTQCD_, GeV, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<MEqq2gZ2ffPowhegQED,Energy> interfaceMinimumpTQEDII
    ("MinimumpTQEDII",
     "The minimum pT for the hard QED emission (initial-initial)",
     &MEqq2gZ2ffPowhegQED::minpTQEDII_, GeV, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<MEqq2gZ2ffPowhegQED,Energy> interfaceMinimumpTQEDIF
    ("MinimumpTQEDIF",
     "The minimum pT for the hard QED emission (initial-final)",
     &MEqq2gZ2ffPowhegQED::minpTQEDIF_, GeV, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<MEqq2gZ2ffPowhegQED,Energy> interfaceMinimumpTQEDFF
    ("MinimumpTQEDFF",
     "The minimum pT for the hard QED emission (final-final)",
     &MEqq2gZ2ffPowhegQED::minpTQEDFF_, GeV, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Reference<MEqq2gZ2ffPowhegQED,ShowerAlpha> interfaceAlphaQCD
    ("AlphaQCD",
     "The object for the calculation of the strong coupling"
     " for the hardest emission",
     &MEqq2gZ2ffPowhegQED::alphaQCD_, false, false, true, false, false);

  static Reference<MEqq2gZ2ffPowhegQED,ShowerAlpha> interfaceAlphaQED
    ("AlphaQED",
     "The object for the calculation of the electromagnetic coupling"
     " for the hardest emission",
     &MEqq2gZ2ffPowhegQED::alphaQED_, false, false, true, false, false);

  static Switch<MEqq2gZ2ffPowhegQED,int> interfaceProcess
    ("Process",
     "Which process to included",
     &MEqq2gZ2ffPowhegQED::process_, 2, false, false);
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
     "Nu_mu",
     "Only include nu_mu nu_mubar as outgoing particles",
     9);
  static SwitchOption interfaceProcessnu_tau
    (interfaceProcess,
     "Nu_tau",
     "Only include nu_tau nu_taubar as outgoing particles",
     10);

  static Parameter<MEqq2gZ2ffPowhegQED,int> interfaceMaxFlavour
    ("MaxFlavour",
     "The maximum flavour of the incoming quarks",
     &MEqq2gZ2ffPowhegQED::maxFlavour_, 5, 1, 5,
     false, false, Interface::limited);

  static Parameter<MEqq2gZ2ffPowhegQED,double> interfacezPower
    ("zPower",
     "The sampling power for z",
     &MEqq2gZ2ffPowhegQED::zPow_, 0.5, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<MEqq2gZ2ffPowhegQED,double> interfaceyPower
    ("yPower",
     "The sampling power for y",
     &MEqq2gZ2ffPowhegQED::yPow_, 0.9, 0.0, 1.0,
     false, false, Interface::limited);

  static Switch<MEqq2gZ2ffPowhegQED,unsigned int> interfaceQEDContributions
    ("QEDContributions",
     "Which QED corrections to include",
     &MEqq2gZ2ffPowhegQED::QEDContributions_, 0, false, false);
  static SwitchOption interfaceQEDContributionsAll
    (interfaceQEDContributions,
     "All",
     "Include all the corrections",
     0);
  static SwitchOption interfaceQEDContributionsISR
    (interfaceQEDContributions,
     "ISR",
     "Only include initial-state QED radiation",
     1);
  static SwitchOption interfaceQEDContributionsFSR
    (interfaceQEDContributions,
     "FSR",
     "Only include final-state QED radiation",
     2);
  static SwitchOption interfaceQEDContributionsISRFSR
    (interfaceQEDContributions,
     "ISRPlusFSR",
     "Only include inital- and final-state QED radiation, no intereference",
     3);

  static Switch<MEqq2gZ2ffPowhegQED,unsigned int> interfaceDipoleSum
    ("DipoleSum",
     "Organization of QED dipoles subtraction",
     &MEqq2gZ2ffPowhegQED::DipoleSum_, 0, false, false);

  static SwitchOption interfaceDipoleSum0
    (interfaceDipoleSum,
     "0",
     "(R-Sum(D-))*|D+|/Sum|D|",
     0);
  static SwitchOption interfaceDipoleSum1
    (interfaceDipoleSum,
     "1",
     "R*|D|/Sum|D|",
     1);

  static Reference<MEqq2gZ2ffPowhegQED,OneLoopFFAWZVertex> interfaceOneLoopVertex
    ("OneLoopVertex",
     "The vertex include the one-loop corrections",
     &MEqq2gZ2ffPowhegQED::oneLoopVertex_, false, false, true, false, false);

  static Switch<MEqq2gZ2ffPowhegQED,int> interfaceQCDRadiationType
    ("QCDRadiationType",
     "Whic types of QCD radiation ot include",
     &MEqq2gZ2ffPowhegQED::QCDRadiationType_, -1, false, false);
  static SwitchOption interfaceQCDRadiationTypeAll
    (interfaceQCDRadiationType,
     "All",
     "All type",
     -1);
  static SwitchOption interfaceQCDRadiationTypeQuarkAntiQuark
    (interfaceQCDRadiationType,
     "QuarkAntiQuark",
     "quark radiates gluon",
     0);
  static SwitchOption interfaceQCDRadiationTypeAntiQuarkQuark
    (interfaceQCDRadiationType,
     "AntiQuarkQuark",
     "antiquark radaites gluon",
     1);
  static SwitchOption interfaceQCDRadiationTypeGluonAntiquark
    (interfaceQCDRadiationType,
     "GluonAntiquark",
     "gluon splits radiation antiquark",
     2);
  static SwitchOption interfaceQCDRadiationTypeQuarkGluon
    (interfaceQCDRadiationType,
     "QuarkGluon",
     "Gluon splits radiating quark",
     3);

}

double MEqq2gZ2ffPowhegQED::subtractedVirtual() const {
  double output(0.);
  // ISR QCD correction
  if(corrections_==1||corrections_==4||corrections_==5)
    output += CFfact_*2.;
  // ISR QED correction
  if((corrections_==2||corrections_==4||corrections_==5) && QEDContributions_ != 2)
    output += sqr(double(mePartonData()[0]->iCharge())/3.)*EMfact_*2.;
  // FSR QED correction
  if((corrections_==2||corrections_==4||corrections_==5) && QEDContributions_ != 1 &&
     mePartonData()[2]->charged()) {
    double mu2 = sqr(mePartonData()[2]->mass())/sHat();
    double mu = sqrt(mu2), mu4 = sqr(mu2), lmu = log(mu);
    double v = sqrt(1.-4.*mu2),v2 = sqr(v),omv = 4.*mu2/(1.+v);
    double f1,f2,fNS,VNS;
    double r = omv/(1.+v),lr(log(r));
    // normal form
    if(mu>1e-4) {
      f1 = 
	( +1. + 3.*log(0.5*(1.+v)) - 1.5*log(0.5*(1.+v2)) + sqr(Constants::pi)/6.
	  - 0.5*sqr(lr) - (1.+v2)/v*(lr*log(1.+v2) + sqr(Constants::pi)/12. 
				     -0.5*log(4.*mu2)*lr + 0.25*sqr(lr)));
      fNS =  -0.5*(1.+2.*v2)*lr/v + 1.5*lr - 2./3.*sqr(Constants::pi) + 0.5*sqr(lr)
	+ (1.+v2)/v*(Herwig::Math::ReLi2(r) + sqr(Constants::pi)/3. -
		     0.25*sqr(lr) + lr*log((2.*v/ (1.+v))));
      VNS = 1.5*log(0.5*(1.+v2)) 
	+ 0.5*(1.+v2)/v*( 2.*lr*log(2.*(1.+v2)/sqr(1.+v))  + 2.*Herwig::Math::ReLi2(sqr(r)) 
			  - 2.*Herwig::Math::ReLi2(2.*v/(1.+v)) - sqr(Constants::pi)/6.)
	+ log(1.-mu) - 2.*log(1.-2.*mu) - 4.*mu2/(1.+v2)*log(mu/(1.-mu)) - mu/(1.-mu)
	+ 4.*(2.*mu2-mu)/(1.+v2) + 0.5*sqr(Constants::pi); 
      f2 = mu2*lr/v;
    }
    // small mass limit
    else {
      f1 = -1./6.*
	( - 6. - 24.*lmu*mu2 - 15.*mu4 - 12.*mu4*lmu - 24.*mu4*sqr(lmu) 
	  + 2.*mu4*sqr(Constants::pi) - 12.*mu2*mu4 - 96.*mu2*mu4*sqr(lmu) 
	  + 8.*mu2*mu4*sqr(Constants::pi) - 80.*mu2*mu4*lmu);
      fNS = - mu2/18.*( + 36.*lmu - 36. - 45.*mu2 + 216.*lmu*mu2 - 24.*mu2*sqr(Constants::pi) 
			+ 72.*mu2*sqr(lmu) - 22.*mu4 + 1032.*mu4 * lmu
			- 96.*mu4*sqr(Constants::pi) + 288.*mu4*sqr(lmu));
      VNS = - mu2/1260.*(-6930. + 7560.*lmu + 2520.*mu - 16695.*mu2 + 1260.*mu2*sqr(Constants::pi) 
			 + 12600.*lmu*mu2 + 1344.*mu*mu2 - 52780.*mu4 + 36960.*mu4*lmu 
			 + 5040.*mu4*sqr(Constants::pi) - 12216.*mu*mu4);
      f2 = mu2*( 2.*lmu + 4.*mu2*lmu + 2.*mu2 + 12.*mu4*lmu + 7.*mu4);
    }
    // subtracted virtual correction
    output += sqr(double(mePartonData()[2]->iCharge())/3.)*
      2.*EMfact_*(f1+fNS+VNS + f2*f2term_/loME_);
  }
  // ISR/FSR interference
  if((corrections_==2||corrections_==4||corrections_==5) && QEDContributions_==0 &&
     mePartonData()[2]->charged()) {
    output += alphaEM_*oneLoopVertex_->
      neutralCurrentQEDME(mePartonData()[0],mePartonData()[1],
       			  mePartonData()[2],mePartonData()[3],
      			  sHat(),tHat(),uHat());
  }
  // genuine EW corrections
  if(corrections_==3||corrections_==5) output += EWterm_;
  // return the total
  return output;
}

double MEqq2gZ2ffPowhegQED::NLOWeight() const {
  // if leading-order return 1
  if(contrib_==0) {
    weights_.resize(1,loME_);
    return loME_;
  }
  // strong coupling
  Energy2 mu2(scale());
  if(!fixedCouplings_) {
    alphaS_  = SM().alphaS (mu2);
    alphaEM_ = SM().alphaEM(mu2);
  }
  // prefactors
  CFfact_ = 4./3.*alphaS_ /Constants::twopi;
  TRfact_ = 1./2.*alphaS_ /Constants::twopi;
  EMfact_ =       alphaEM_/Constants::twopi;
  // virtual pieces
  double virt = subtractedVirtual();
  // extract the partons and stuff for the real emission
  //  and collinear counter terms
  // hadrons
  pair<tcBeamPtr,tcBeamPtr> hadrons= 
    make_pair(dynamic_ptr_cast<tcBeamPtr>(lastParticles().first->dataPtr() ),
	      dynamic_ptr_cast<tcBeamPtr>(lastParticles().second->dataPtr()));
  // momentum fractions
  pair<double,double> x = make_pair(lastX1(),lastX2());
  // partons
  pair<tcPDPtr,tcPDPtr> partons = make_pair(mePartonData()[0],mePartonData()[1]);
  // If necessary swap the particle data objects so that 
  // first beam gives the incoming quark
  if(lastPartons().first ->dataPtr()!=partons.first) {
    swap(x.first,x.second);
    swap(hadrons.first,hadrons.second);
  }
  // convert the values of z tilde to z
  pair<double,double> z;
  pair<double,double> zJac;
  double rhomax(pow(1.-x.first,1.-power_));
  double rho = zTilde_*rhomax;
  z.first = 1.-pow(rho,1./(1.-power_));
  zJac.first = rhomax*pow(1.-z.first,power_)/(1.-power_);
  rhomax = pow(1.-x.second,1.-power_);
  rho = zTilde_*rhomax; 
  z.second = 1.-pow(rho,1./(1.-power_));
  zJac.second = rhomax*pow(1.-z.second,power_)/(1.-power_);
  // calculate the PDFs
  pair<double,double> oldqPDF = 
    make_pair(hadrons.first ->pdf()->xfx(hadrons.first ,partons.first ,scale(),
					 x.first )/x.first ,
	      hadrons.second->pdf()->xfx(hadrons.second,partons.second,scale(),
					 x.second)/x.second);
  // real/coll q/qbar
  pair<double,double> newqPDF = 
    make_pair(hadrons.first ->pdf()->xfx(hadrons.first ,partons.first ,scale(),
					 x.first /z.first )*z.first /x.first ,
	      hadrons.second->pdf()->xfx(hadrons.second,partons.second,scale(),
					 x.second/z.second)*z.second/x.second);
  // real/coll gluon
  pair<double,double> newgPDF =  
    make_pair(hadrons.first ->pdf()->xfx(hadrons.first ,gluon_,scale(),
					 x.first /z.first )*z.first /x.first ,
	      hadrons.second->pdf()->xfx(hadrons.second,gluon_,scale(),
					 x.second/z.second)*z.second/x.second);
  pair<double,double> newpPDF(make_pair(-1.,-1.));
  // real/coll photon
  if(incomingPhotons_) {
    // if check the PDF has photons
    cPDVector partons = hadrons.first ->pdf()->partons(hadrons.first);
    if(find(partons.begin(),partons.end(),gamma_)!=partons.end())
      newpPDF.first  = hadrons.first ->pdf()->xfx(hadrons.first ,gamma_,scale(),
						  x.first /z.first )*z.first /x.first;
    partons = hadrons.second->pdf()->partons(hadrons.second);
    if(find(partons.begin(),partons.end(),gamma_)!=partons.end())
      newpPDF.second = hadrons.second->pdf()->xfx(hadrons.second,gluon_,scale(),
						  x.second/z.second)*z.second/x.second;
  }
  // collinear remnants 
  double coll = 0.;
  // common terms
  // q -> q
  double collQQ       = collinearQuark(x.first ,mu2,zJac.first ,z.first ,
				       oldqPDF.first ,newqPDF.first );
  // qbar -> qbar
  double collQbarQbar = collinearQuark(x.second,mu2,zJac.second,z.second,
				       oldqPDF.second,newqPDF.second);
  // QCD
  if(corrections_==1||corrections_==4||corrections_==5) {
    // g -> q
    double collGQ       = collinearBoson(mu2,zJac.first ,z.first ,
					 oldqPDF.first ,newgPDF.first );
    // g -> qbar
    double collGQbar    = collinearBoson(mu2,zJac.second,z.second,
					 oldqPDF.second,newgPDF.second);
    // add to sum
    coll += CFfact_*(collQQ+collQbarQbar)+TRfact_*(collGQ+collGQbar);
  }
  // QED
  if((corrections_==2||corrections_==4||corrections_==5) &&
     QEDContributions_!=2 ) {
    // gamma -> q
    double collPQ       = newpPDF.first  < 0. ? 0. : 
      collinearBoson(mu2,zJac.first ,z.first ,oldqPDF.first ,newpPDF.first );
    // gamma -> qbar
    double collPQbar    = newpPDF.second < 0. ? 0. : 
      collinearBoson(mu2,zJac.second,z.second,oldqPDF.second,newpPDF.second);
    // add to sum
    coll += sqr(double(mePartonData()[0]->iCharge())/3.)*EMfact_*
      (collQQ+collQbarQbar+3.*(collPQ+collPQbar));

    // q -> q (DIS factorization scheme)
    double collQQDIS       = collinearQuarkKDIS(x.first,zJac.first ,z.first ,
						oldqPDF.first ,newqPDF.first );
    // qbar -> qbar (DIS factorization scheme)
    double collQbarQbarDIS = collinearQuarkKDIS(x.second,zJac.second,z.second,
						oldqPDF.second,newqPDF.second );
    // gamma -> q (DIS factorization scheme)
    double collPQDIS       = newpPDF.first  < 0. ? 0. : 
      collinearBosonKDIS(zJac.first ,z.first ,oldqPDF.first ,newpPDF.first );
    // gamma -> qbar (DIS factorization scheme)
    double collPQbarDIS    = newpPDF.second < 0. ? 0. : 
      collinearBosonKDIS(zJac.second,z.second,oldqPDF.second,newpPDF.second);
    // add to sum
    coll+= sqr(double(mePartonData()[0]->iCharge())/3.)*EMfact_*
      (collQQDIS+collQbarQbarDIS+3.*(collPQDIS+collPQbarDIS));

    // q -> q (QED IF piece)
    double collQIF       = collinearQuarkIF(x.first,zJac.first ,z.first ,
					    oldqPDF.first ,newqPDF.first ) ;
    // qbar -> qbar (QED IF piece)
    double collQbarIF    = collinearQuarkIF(x.second,zJac.second,z.second,
					    oldqPDF.second,newqPDF.second) ;
    // SumPop_IF corresponds to CS, eq. 10.25 (Ti*Tap terms of the P operator), 
    // for the QED case.
    // Notice that for Z production the analogous terms in the K operator (CS, eq. 10.24)
    // vanish when the sum is performed.
    double SumPop_IF=0.;
    for (unsigned int iin=0; iin<2; ++iin) {
      double collIF = iin==0 ? collQIF : collQbarIF ;
      for(unsigned int iout=2; iout<4; ++iout) {
	// - sign because outgoing particles should be counted with a minus //ER
	int QiQo= - mePartonData()[iin]->iCharge()*mePartonData()[iout]->iCharge(); 
	if( QiQo != 0) {  
	  Energy2 sioHat = 2.*meMomenta()[iin].dot(meMomenta()[iout]) ;
	  //cout<<iin<<iout<<"  "<<collIF/collQIF<<" "<<collIF/collQbarIF<<endl;
	  //Energy2 sioHat2 = ((iin+iout) %2 ==0 ) ? -tHat() : -uHat(); 
	  //cout << sioHat/sioHat2 << endl;
	  SumPop_IF += double(QiQo)/9. * log(mu2/sioHat) * collIF;
	}
      }
    }
    // add to sum
    coll += SumPop_IF*EMfact_ ;
  }
  // add up the virtual and remnant terms
  double wgt = loME_*( 1. + virt + coll );
  // real QCD emission
  vector<double> realQCD1,realQCD2;
  if(corrections_==1||corrections_==4||corrections_==5) {
    realQCD1 = subtractedRealQCD(x,z. first,zJac. first,
				 oldqPDF. first,newqPDF. first,
				 newgPDF. first, II12);
    realQCD2 = subtractedRealQCD(x,z.second,zJac.second, 
				 oldqPDF.second,newqPDF.second,
				 newgPDF.second, II21);
    wgt += realQCD1[0] + realQCD1[2] +realQCD2[0] +realQCD2[2];
  }
  // real QED emission
  vector<double> realQED1(4,0.),realQED2(4,0.),
    realQED3(2,0.),realQED4(2,0.),
    realQED5(2,0.),realQED6(2,0.),
    realQED7(2,0.),realQED8(2,0.);
  if(corrections_==2||corrections_==4||corrections_==5) {
    // ISR
    if(QEDContributions_!=2) {
//       cout<<"====II======="<<endl;
//       cout<<mePartonData()[0]->PDGName()<<" "<<
// 	mePartonData()[1]->PDGName()<<" "<<
// 	mePartonData()[2]->PDGName()<<" "<<
// 	mePartonData()[3]->PDGName()<<" "<<endl;

      realQED1 = subtractedRealQED(x,z. first,zJac. first,
				   oldqPDF. first,newqPDF. first,
				   newpPDF. first, II12);
      realQED2 = subtractedRealQED(x,z.second,zJac.second, 
				   oldqPDF.second,newqPDF.second,
				   newpPDF.second, II21);
      wgt += realQED1[0] + realQED1[2] +realQED2[0] +realQED2[2];
    }
    // FSR
    if(QEDContributions_!=1) {
//       cout<<"====FF======="<<endl;
//       cout<<mePartonData()[0]->PDGName()<<" "<<
// 	mePartonData()[1]->PDGName()<<" "<<
// 	mePartonData()[2]->PDGName()<<" "<<
// 	mePartonData()[3]->PDGName()<<" "<<endl;

      realQED3 = subtractedRealQED(x,zTilde_,1.,oldqPDF. first,
				   oldqPDF. first,0.,FF34);
      realQED4 = subtractedRealQED(x,zTilde_,1.,oldqPDF. first,
				   oldqPDF. first,0.,FF43);
      wgt += realQED3[0] + realQED4[0];
    }
    // ISR/FSR interference
    if(QEDContributions_==0) {
//       cout<<"====IF======="<<endl;
//       cout<<mePartonData()[0]->PDGName()<<" "<<
// 	mePartonData()[1]->PDGName()<<" "<<
// 	mePartonData()[2]->PDGName()<<" "<<
// 	mePartonData()[3]->PDGName()<<" "<<endl;
      
      if(DipoleSum_==0 ||
	 mePartonData()[0]->iCharge()*mePartonData()[2]->iCharge()>0)
	realQED5 = subtractedRealQED(x,z. first,zJac. first,
				     oldqPDF. first,newqPDF. first,
				     newpPDF. first, IF13);
      if(DipoleSum_==0 ||
	 mePartonData()[0]->iCharge()*mePartonData()[3]->iCharge()>0)
	realQED6 = subtractedRealQED(x,z. first,zJac. first,
				     oldqPDF. first,newqPDF. first,
				     newpPDF. first, IF14);
      if(DipoleSum_==0 ||
	 mePartonData()[1]->iCharge()*mePartonData()[2]->iCharge()>0)
	realQED7 = subtractedRealQED(x,z.second,zJac.second, 
				     oldqPDF.second,newqPDF.second,
				     newpPDF.second, IF23);
      if(DipoleSum_==0 ||
	 mePartonData()[1]->iCharge()*mePartonData()[3]->iCharge()>0)
	realQED8 = subtractedRealQED(x,z.second,zJac.second, 
				     oldqPDF.second,newqPDF.second,
				     newpPDF.second, IF24);
      //       cout<<"realqed= "<<realQED5[0]<<" "<<realQED6[0]<<" "<<realQED7[0]<<" "<<realQED8[0]<<" "<<endl;
      wgt += (realQED5[0] + realQED6[0] + realQED7[0] + realQED8[0] );
    }
  }
  if(isnan(wgt)||isinf(wgt)) {
    generator()->log() << "testing bad weight "
		       << virt << " " << coll << "\n";
    generator()->log() << "QCD1 ";
    for(unsigned int ix=0;ix<realQCD1.size();++ix)
      generator()->log() << realQCD1[ix] << " ";
    generator()->log() << "\n";
    generator()->log() << "QCD2 ";
    for(unsigned int ix=0;ix<realQCD2.size();++ix)
      generator()->log() << realQCD2[ix] << " ";
    generator()->log() << "\n";


    generator()->log() << "QED1 ";
    for(unsigned int ix=0;ix<realQED1.size();++ix)
      generator()->log() << realQED1[ix] << " ";
    generator()->log() << "\n";
    generator()->log() << "QED2 ";
    for(unsigned int ix=0;ix<realQED2.size();++ix)
      generator()->log() << realQED2[ix] << " ";
    generator()->log() << "\n";
    generator()->log() << "QED3 ";
    for(unsigned int ix=0;ix<realQED3.size();++ix)
      generator()->log() << realQED3[ix] << " ";
    generator()->log() << "\n";
    generator()->log() << "QED4 ";
    for(unsigned int ix=0;ix<realQED4.size();++ix)
      generator()->log() << realQED4[ix] << " ";
    generator()->log() << "\n";
    generator()->log() << "QED5 ";
    for(unsigned int ix=0;ix<realQED5.size();++ix)
      generator()->log() << realQED5[ix] << " ";
    generator()->log() << "\n";
    generator()->log() << "QED6 ";
    for(unsigned int ix=0;ix<realQED6.size();++ix)
      generator()->log() << realQED6[ix] << " ";
    generator()->log() << "\n";
    generator()->log() << "QED7 ";
    for(unsigned int ix=0;ix<realQED7.size();++ix)
      generator()->log() << realQED7[ix] << " ";
    generator()->log() << "\n";
    generator()->log() << "QED8 ";
    for(unsigned int ix=0;ix<realQED8.size();++ix)
      generator()->log() << realQED8[ix] << " ";
    generator()->log() << "\n";
    generator()->log() << "testing z " << z.first << " " << z.second << "\n";
    generator()->log() << "testing z " << 1.-z.first << " " << 1.-z.second << "\n";
    generator()->log() << flush;
    assert(false);
  }
  weights_.resize(15,0.);
  if(corrections_==1||corrections_==4||corrections_==5) {
    weights_[ 1] = realQCD1[1];
    weights_[ 2] = realQCD1[3];
    weights_[ 3] = realQCD2[1];
    weights_[ 4] = realQCD2[3];
  }
  if(corrections_==2||corrections_==4||corrections_==5) {
    weights_[ 5] = realQED1[1];
    weights_[ 6] = realQED1[3];
    weights_[ 7] = realQED2[1];
    weights_[ 8] = realQED2[3];
    weights_[ 9] = realQED3[1];
    weights_[10] = realQED4[1];
    weights_[11] = realQED5[1];
    weights_[12] = realQED6[1];
    weights_[13] = realQED7[1];
    weights_[14] = realQED8[1];
  }
  if( contrib_ < 3 ) {
    weights_[0] = wgt;
    wgt = std::accumulate(weights_.begin(),weights_.end(),0.);
    return contrib_ == 1 ? max(0.,wgt) : max(0.,-wgt);
  }
  else if(contrib_==3) {
    weights_[0] = 0.;
    return std::accumulate(weights_.begin(),weights_.end(),0.);
  }
  else 
    return wgt;
}

double 
MEqq2gZ2ffPowhegQED::collinearQuark(double x, Energy2 mu2, 
				    double jac, double z,
				    double oldPDF, double newPDF) const {
  if(1.-z < 1.e-8) return 0.;
  return 
    // this bit is multiplied by LO PDF
    sqr(Constants::pi)/3.-5.+2.*sqr(log(1.-x ))
    +(1.5+2.*log(1.-x ))*log(sHat()/mu2)
    // NLO PDF bit
    +jac /z * newPDF /oldPDF *
    (1.-z -(1.+z )*log(sqr(1.-z )/z )
     -(1.+z )*log(sHat()/mu2)-2.*log(z )/(1.-z ))
    // + function bit
    +jac /z *(newPDF /oldPDF -z )*
    2./(1.-z )*log(sHat()*sqr(1.-z )/mu2);
}

double 
MEqq2gZ2ffPowhegQED::collinearBoson(Energy2 mu2, double jac, double z,
				    double oldPDF, double newPDF) const {
  if(1.-z < 1.e-8) return 0.;
  return jac/z*newPDF/oldPDF*
    ((sqr(z)+sqr(1.-z))*log(sqr(1.-z)*sHat()/z/mu2)
     +2.*z*(1.-z));
}

double 
MEqq2gZ2ffPowhegQED::collinearQuarkIF(double x, double jac, double z,
				    double oldPDF, double newPDF) const {
  if(1.-z < 1.e-8) return 0.;
  return 
    // this bit is multiplied by LO PDF
    1.5+2*log(1.-x)
    // plus distribution bit
    +jac /(1.-z) * 
    (newPDF /oldPDF *(1+z*z) /z - 2) ;
}

double 
MEqq2gZ2ffPowhegQED::collinearQuarkKDIS(double x, double jac, double z,
				    double oldPDF, double newPDF) const {
  // this is the function defined in eq. C.25 of CS
  double fun=(1+z*z)/(1-z)*(log((1.-z)/z)-0.75) + (9+5*z)/4.;
  // this is artanh(1-2x)
  double atanhom2x=0.5*log((1.-x)/x);
  // this is \int_0^x (fun(z) dz)
  double intfun0x=0.5*(4*Math::ReLi2(x)
		       + x *(2*x+7) 
		       -2 *log(1.-x) *(log(1.-x) -2*log(x) -3)
		       -2 *x*(x+2)*atanhom2x) ;
  if(1.-z < 1.e-8) return 0.;
  return 
    // this bit is multiplied by LO PDF
    intfun0x
    // plus distribution bit
    - jac * fun *
    (newPDF /oldPDF /z -1);
}

double 
MEqq2gZ2ffPowhegQED::collinearBosonKDIS(double jac, double z,
				    double oldPDF, double newPDF) const {
  if(1.-z < 1.e-8) return 0.;
  return jac/z*newPDF/oldPDF*
    ((sqr(z)+sqr(1.-z))*log((1.-z)/z) + 8*z*(1.-z) -1.);
}




vector<double> 
MEqq2gZ2ffPowhegQED::subtractedRealQCD(pair<double,double> x, double z,
				       double zJac, double oldqPDF,
				       double newqPDF, double newgPDF,
				       DipoleType dipole) const {
  double vt   = vTilde_*(1.-z);
  double vJac = 1.-z;
  Energy pT   = sqrt(sHat()*vt*(1.-vt-z)/z);
  // rapidities
  double rapidity;
  if(dipole==II12) {
    rapidity = -log(x.second*sqrt(lastS())/pT*vt);
  }
  else if(dipole==II21) {
    rapidity =  log(x.first *sqrt(lastS())/pT*vt);
  }
  else {
    assert(false);
  }
  // CMS system
  Energy rs=sqrt(lastS());
  Lorentz5Momentum pcmf = Lorentz5Momentum(ZERO,ZERO,0.5*rs*(x.first-x.second),
					   0.5*rs*(x.first+x.second));
  pcmf.rescaleMass();
  Boost blab(pcmf.boostVector());
  // emission from the quark radiation
  vector<Lorentz5Momentum> pnew(5);
  if(dipole==II12) {
    pnew [0] = Lorentz5Momentum(ZERO,ZERO,0.5*rs*x.first/z,
				0.5*rs*x.first/z,ZERO);
    pnew [1] = Lorentz5Momentum(ZERO,ZERO,-0.5*rs*x.second,
				0.5*rs*x.second,ZERO) ;
  }
  else if(dipole==II21) {
    pnew[0] = Lorentz5Momentum(ZERO,ZERO,0.5*rs*x.first,
			       0.5*rs*x.first,ZERO);
    pnew[1] = Lorentz5Momentum(ZERO,ZERO,-0.5*rs*x.second/z,
			       0.5*rs*x.second/z,ZERO) ;
  }
  else {
    assert(false);
  }
  pnew [2] = meMomenta()[2];
  pnew [3] = meMomenta()[3];
  pnew [4] = Lorentz5Momentum(pT*cos(phi_),pT*sin(phi_),
			      pT*sinh(rapidity),
			      pT*cosh(rapidity), ZERO);
  Lorentz5Momentum K  = pnew [0]+pnew [1]-pnew [4];
  Lorentz5Momentum Kt = pcmf;
  Lorentz5Momentum Ksum = K+Kt;
  Energy2 K2 = K.m2();
  Energy2 Ksum2 = Ksum.m2();
  for(unsigned int ix=2;ix<4;++ix) {
    pnew [ix].boost(blab);
    pnew [ix] = pnew [ix] - 2.*Ksum*(Ksum*pnew [ix])/Ksum2
      +2*K*(Kt*pnew [ix])/K2;
  }
  // phase-space prefactors
  // double phase = zJac*vJac/sqr(z);
  double phase = zJac*vJac/z;
  // real emission q qbar
  vector<double> output(4,0.);
  if(dipole==II12) {
    realEmissionQCDGluon1_ = pnew;
    realEmissionQCDQuark1_ = pnew;
  }
  else if(dipole==II21) {
    realEmissionQCDGluon2_ = pnew;
    realEmissionQCDQuark2_ = pnew;
  }
  else {
    assert(false);
  }
  if(!(zTilde_<1e-7 || vt<1e-7 || 1.-z-vt < 1e-7 )) {
    pair<double,double> realQQ = subtractedQCDMEqqbar(pnew,dipole,true);
    double fact1 = CFfact_*phase*newqPDF/oldqPDF;
    pair<double,double> realGQ = subtractedQCDMEgqbar(pnew,dipole,true);
    double fact2 = TRfact_*phase*newgPDF/oldqPDF;
    output[0] = realQQ.first *fact1;
    output[1] = realQQ.second*fact1;
    output[2] = realGQ.first *fact2;
    output[3] = realGQ.second*fact2;
  }
  // return the answer
  return output;
}

pair<double,double> 
MEqq2gZ2ffPowhegQED::subtractedQCDMEqqbar(const vector<Lorentz5Momentum> & p,
					  DipoleType dipole,bool subtract) const {
  // use the inheriting class to calculate the matrix element
  cPDVector particles(mePartonData());
  particles.push_back(gluon_);
  double me = 0.75*realQCDME(particles,p);
  // compute the two dipole terms
  double x = (p[0]*p[1]-p[4]*p[1]-p[4]*p[0])/(p[0]*p[1]);
  Lorentz5Momentum Kt = p[0]+p[1]-p[4];
  vector<Lorentz5Momentum> pa(4),pb(4);
  // momenta for q -> q g emission
  pa[0] = x*p[0];
  pa[1] =   p[1];
  Lorentz5Momentum K = pa[0]+pa[1];
  Lorentz5Momentum Ksum = K+Kt;
  Energy2 K2 = K.m2();
  Energy2 Ksum2 = Ksum.m2();
  for(unsigned int ix=2;ix<4;++ix)
    pa[ix] = p[ix]-2.*Ksum*(Ksum*p[ix])/Ksum2+2*K*(Kt*p[ix])/K2;
  // first LO matrix element
  double lo1 = loME(mePartonData(),pa,false);
  // momenta for qbar -> qbar g emission
  pb[0] =   p[0];
  pb[1] = x*p[1];
  K = pb[0]+pb[1];
  Ksum = K+Kt;
  K2 = K.m2();
  Ksum2 = Ksum.m2();
  for(unsigned int ix=2;ix<4;++ix)
    pb[ix] = p[ix]-2.*Ksum*(Ksum*p[ix])/Ksum2+2*K*(Kt*p[ix])/K2;
  // second LO matrix element
  double lo2 = loME(mePartonData(),pb,false);
  // first dipole
  InvEnergy2 D1 = 0.5/(p[0]*p[4])/x*(2./(1.-x)-(1.+x));
  // second dipole
  InvEnergy2 D2 = 0.5/(p[1]*p[4])/x*(2./(1.-x)-(1.+x));
  // results
  pair<double,double> supressionFactor = supressionFunction(sqr(p[4].x())+sqr(p[4].y()));
  pair<double,double> output = make_pair(0.,0.);
  if(lo1>0.&&lo2>0.) {
    if(dipole==II12) {
      me *= abs(D1)*lo1/(abs(D1)*lo1+abs(D2)*lo2);
      if(subtract) {
	output.first  = sHat()*(UnitRemoval::InvE2*me*supressionFactor.first -D1*lo1);
      }
      else {
	output.first  = sHat()*UnitRemoval::InvE2*me*supressionFactor.first;
      }
      output.second = sHat()*(UnitRemoval::InvE2*me*supressionFactor.second);
    }
    else if(dipole==II21) {
      me *= abs(D2)*lo2/(abs(D1)*lo1+abs(D2)*lo2);
      if(subtract) {
	output.first  = sHat()*(UnitRemoval::InvE2*me*supressionFactor.first -D2*lo2);
      }
      else {
	output.first  = sHat()*UnitRemoval::InvE2*me*supressionFactor.first;
      }
      output.second = sHat()*(UnitRemoval::InvE2*me*supressionFactor.second);
    }
    else {
      assert(false);
    }
  }
  return output;
}

pair<double,double>
MEqq2gZ2ffPowhegQED::subtractedQCDMEgqbar(const vector<Lorentz5Momentum> & p,
					  DipoleType dipole,bool subtract) const {
  // use the inheriting class to calculate the matrix element
  cPDVector particles(mePartonData());
  if(dipole==II12) {
    particles.push_back(particles[0]->CC());
    particles[0] = gluon_;
  }
  else if(dipole==II21) {
    particles.push_back(particles[1]->CC());
    particles[1] = gluon_;
  }
  else {
    assert(false);
  }
  double me = 2.*realQCDME(particles,p);
  // compute the two dipole terms
  double x = 1.-(p[4]*p[1]+p[4]*p[0])/(p[0]*p[1]);
  Lorentz5Momentum Kt = p[0]+p[1]-p[4];
  vector<Lorentz5Momentum> pa(4);
  // momenta for ISR
  if(dipole==II12) {
    pa[0] = x*p[0];
    pa[1] =   p[1];
  }
  else if(dipole==II21) {
    pa[0] =   p[0];
    pa[1] = x*p[1];
  }
  else {
    assert(false);
  }
  Lorentz5Momentum K = pa[0]+pa[1];
  Lorentz5Momentum Ksum = K+Kt;
  Energy2 K2 = K.m2();
  Energy2 Ksum2 = Ksum.m2();
  for(unsigned int ix=2;ix<4;++ix)
    pa[ix] = p[ix]-2.*Ksum*(Ksum*p[ix])/Ksum2+2*K*(Kt*p[ix])/K2; 
  // first LO matrix element 
  double lo1 = loME(mePartonData(),pa,false); 
  // dipole
  InvEnergy2 D1;
  if(dipole==II12) {
    D1 =  0.5/(p[0]*p[4])/x*(1.-2.*x*(1.-x));
  }
  else if(dipole==II21) {
    D1 =  0.5/(p[1]*p[4])/x*(1.-2.*x*(1.-x));
  }
  else {
    assert(false);
  }
  pair<double,double> supressionFactor = 
    supressionFunction(sqr(p[4].x())+sqr(p[4].y()));
  if(subtract) {
    return make_pair(sHat()*(UnitRemoval::InvE2*me*supressionFactor.first -D1*lo1),
		     sHat()*(UnitRemoval::InvE2*me*supressionFactor.second));
  }
  else {
    return make_pair(sHat()*(UnitRemoval::InvE2*me*supressionFactor.first ),
		     sHat()*(UnitRemoval::InvE2*me*supressionFactor.second));
  }
}

double MEqq2gZ2ffPowhegQED::realQCDME(const cPDVector & particles,
				      const vector<Lorentz5Momentum> & momenta) const {
  vector<SpinorWaveFunction>    fin(2),aout(2);
  vector<SpinorBarWaveFunction> ain(2),fout(2);
  vector<VectorWaveFunction> gluon(2);
  // wavefunctions for the q qbar -> l+ l- g process
  if(particles[0]->id()==-particles[1]->id()) {
    for( unsigned int i = 0; i < 2; ++i ) {
      fin[i]   = SpinorWaveFunction   (momenta[0],particles[0],  i,incoming);
      ain[i]   = SpinorBarWaveFunction(momenta[1],particles[1],  i,incoming);
      gluon[i] = VectorWaveFunction   (momenta[4],particles[4],2*i,outgoing);
    }
  }
  else if(particles[0]->id()==ParticleID::g && particles[1]->id()<0) {
    for( unsigned int i = 0; i < 2; ++i ) {
      fin[i]   = SpinorWaveFunction   (momenta[4],particles[4],  i,outgoing);
      ain[i]   = SpinorBarWaveFunction(momenta[1],particles[1],  i,incoming);
      gluon[i] = VectorWaveFunction   (momenta[0],particles[0],2*i,incoming);
    }
  }
  else if(particles[0]->id()>0 && particles[1]->id()==ParticleID::g) {
    for( unsigned int i = 0; i < 2; ++i ) {
      fin[i]   = SpinorWaveFunction   (momenta[0],particles[0],  i,incoming);
      ain[i]   = SpinorBarWaveFunction(momenta[4],particles[4],  i,outgoing);
      gluon[i] = VectorWaveFunction   (momenta[1],particles[1],2*i,incoming);
    }
  }
  else {
    for(unsigned int ix=0;ix<particles.size();++ix) {
      generator()->log() << particles[ix]->PDGName() << " " << momenta[ix]/GeV << "\n";
    }
    assert(false);
  }
  // wavefunctions for the leptons
  for(unsigned int i=0; i<2; ++i) {
    fout[i] = SpinorBarWaveFunction(momenta[2],particles[2],i,outgoing);
    aout[i] = SpinorWaveFunction   (momenta[3],particles[3],i,outgoing);
  }
  double output(0.);
  vector<Complex> diag(4,0.);
  Energy2 q2 = scale();
  for(unsigned int lhel1=0;lhel1<2;++lhel1) {
    for(unsigned int lhel2=0;lhel2<2;++lhel2) {
      VectorWaveFunction Zwave = oneLoopVertex_->evaluate(q2, 1, Z0_,
							  aout[lhel2],fout[lhel1]);
      VectorWaveFunction Pwave = oneLoopVertex_->evaluate(q2, 1, gamma_,
							  aout[lhel2],fout[lhel1]);
      for(unsigned int ihel1=0;ihel1<2;++ihel1) {
	for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	  for(unsigned int ohel1=0;ohel1<2;++ohel1) {
	    //////// Z      diagrams //////////
	    // first Z diagram
	    SpinorWaveFunction inters = FFGVertex_->
	      evaluate(q2,5,fin[ihel1].particle(),fin[ihel1],gluon[ohel1]);
 	    diag[0] = oneLoopVertex_->evaluate(q2, inters, ain[ihel2], Zwave);
	    // second Z diagram
	    SpinorBarWaveFunction interb = FFGVertex_->
	      evaluate(q2,5,ain[ihel1].particle(),ain[ihel2],gluon[ohel1]);
	    diag[1] = oneLoopVertex_->evaluate(q2, fin[ihel1], interb, Zwave);
	    //////// photon diagrams //////////
	    if(particles[2]->id()==-particles[3]->id()&&
	       particles[2]->charged()) {
	      // first photon diagram
	      diag[2] = oneLoopVertex_->evaluate(q2, inters, ain[ihel2], Pwave);
	      // second photon diagram
	      diag[3] = oneLoopVertex_->evaluate(q2, fin[ihel1], interb, Pwave);
	    }
	    // add them up
	    output += norm(std::accumulate(diag.begin(),diag.end(),Complex(0.)));
	  }
	}
      }
    }
  }
  // colour and spin factors
  if(particles[0]->id()==-particles[1]->id()) {
    output *= 1./9.;
  }
  else  {
    output *= 1./24.;
  }
  // divided by 2 g_S^2
  return 0.5*output/norm(FFGVertex_->norm());
}

vector<double> 
MEqq2gZ2ffPowhegQED::subtractedRealQED(pair<double,double> x, double z,
				       double zJac, double oldqPDF,
				       double newqPDF, double newpPDF,
				       DipoleType dipole) const {
  // ISR
  if(dipole==II12||dipole==II21) {
    double vt   = vTilde_*(1.-z);
    double vJac = 1.-z;
    Energy pT   = sqrt(sHat()*vt*(1.-vt-z)/z);
    // rapidities
    double rapidity;
    if(dipole==II12) {
      rapidity = -log(x.second*sqrt(lastS())/pT*vt);
    }
    else if(dipole==II21) {
      rapidity =  log(x.first *sqrt(lastS())/pT*vt);
    }
    else {
      assert(false);
    }
    // CMS system
    Energy rs=sqrt(lastS());
    Lorentz5Momentum pcmf = Lorentz5Momentum(ZERO,ZERO,0.5*rs*(x.first-x.second),
					     0.5*rs*(x.first+x.second));
    pcmf.rescaleMass();
    Boost blab(pcmf.boostVector());
    // emission from the quark radiation
    vector<Lorentz5Momentum> pnew(5);
    if(dipole==II12) {
      pnew [0] = Lorentz5Momentum(ZERO,ZERO,0.5*rs*x.first/z,
				  0.5*rs*x.first/z,ZERO);
      pnew [1] = Lorentz5Momentum(ZERO,ZERO,-0.5*rs*x.second,
				  0.5*rs*x.second,ZERO) ;
    }
    else if(dipole==II21) {
      pnew[0] = Lorentz5Momentum(ZERO,ZERO,0.5*rs*x.first,
				 0.5*rs*x.first,ZERO);
      pnew[1] = Lorentz5Momentum(ZERO,ZERO,-0.5*rs*x.second/z,
				 0.5*rs*x.second/z,ZERO) ;
    }
    else {
      assert(false);
    }
    pnew [2] = meMomenta()[2];
    pnew [3] = meMomenta()[3];
    pnew [4] = Lorentz5Momentum(pT*cos(phi_),pT*sin(phi_),
				pT*sinh(rapidity),
				pT*cosh(rapidity), ZERO);
    Lorentz5Momentum K  = pnew [0]+pnew [1]-pnew [4];
    Lorentz5Momentum Kt = pcmf;
    Lorentz5Momentum Ksum = K+Kt;
    Energy2 K2 = K.m2();
    Energy2 Ksum2 = Ksum.m2();
    for(unsigned int ix=2;ix<4;++ix) {
      pnew [ix].boost(blab);
      pnew [ix] = pnew [ix] - 2.*Ksum*(Ksum*pnew [ix])/Ksum2
	+2*K*(Kt*pnew [ix])/K2;
    }
    // phase-space prefactors
    double phase = zJac*vJac/z;
    if(dipole==II12) {
      realEmissionQEDPhoton1_ = pnew;
      realEmissionQEDQuark1_  = pnew;
    }
    else if(dipole==II21) {
      realEmissionQEDPhoton2_ = pnew;
      realEmissionQEDQuark2_  = pnew;
    }
    else {
      assert(false);
    }
    vector<double> output(4,0.);
    if(!(zTilde_<1e-7 || vt<1e-7 || 1.-z-vt < 1e-7 )) {
      // real emission q qbar
      pair<double,double> realQQ = subtractedQEDMEqqbar(pnew,dipole,true);
      double fact1 = EMfact_*phase*newqPDF/oldqPDF;
      pair<double,double> realPQ(make_pair(0.,0.));
      double fact2 = 0.;
      if(newpPDF>0.) {
	realPQ = subtractedQEDMEpqbar(pnew,dipole,true);
	fact2 = EMfact_*phase*newpPDF/oldqPDF;
      }
      output[0] = realQQ.first *fact1;
      output[1] = realQQ.second*fact1;
      output[2] = realPQ.first *fact2;
      output[3] = realPQ.second*fact2;
    }
    // return the answer
    return output;
  }
  // FSR
  else if(dipole==FF34||dipole==FF43) {
    // reduced mass
    double mu2 = sqr(mePartonData()[2]->mass())/sHat();
    double mu  = sqrt(mu2);
    // jacobian
    double jac = zJac;
    // generate y
    double yminus = 0.; 
    double yplus  = 1.-2.*mu*(1.-mu)/(1.-2*mu2);
    double rhoymax = pow(yplus-yminus,1.-yPow_);
    double rhoy = zTilde_*rhoymax;
    double y = yminus+pow(rhoy,1./(1.-yPow_));
    jac *= pow(y-yminus,yPow_)*rhoymax/(1.-yPow_);
    // generate z 
    double vt = sqrt(max(sqr(2.*mu2+(1.-2.*mu2)*(1.-y))-4.*mu2,0.))/(1.-2.*mu2)/(1.-y);
    double zplus  = (1.+vt)*(1.-2.*mu2)*y/2./(mu2 +(1.-2.*mu2)*y);
    double zminus = (1.-vt)*(1.-2.*mu2)*y/2./(mu2 +(1.-2.*mu2)*y);
    double rhozmax = pow(zplus-zminus,1.-zPow_);
    double rhoz = vTilde_*rhozmax;
    double z = zminus+pow(rhoz,1./(1.-zPow_));
    jac *= pow(z-zminus,zPow_)*rhozmax/(1.-zPow_);
    // calculate x1,x2,x3 and xT 
    double x2 = 1. - y*(1.-2.*mu2);
    double x1 = 1. - z*(x2-2.*mu2);
    double x3 = 2.-x1-x2;
    double xT = sqrt(max(0.,sqr(x3) -0.25*sqr(sqr(x2)+sqr(x3)-sqr(x1))/(sqr(x2)-4.*mu2)));
    // calculate the momenta
    Energy M = sqrt(sHat());
    Lorentz5Momentum pspect(ZERO,ZERO,-0.5*M*sqrt(max(sqr(x2)-4.*mu2,0.)),0.5*M*x2,M*mu); 
    Lorentz5Momentum pemit (-0.5*M*xT*cos(phi_),-0.5*M*xT*sin(phi_),
			    0.5*M*sqrt(max(sqr(x1)-sqr(xT)-4.*mu2,0.)),0.5*M*x1,M*mu);
    Lorentz5Momentum pphoton( 0.5*M*xT*cos(phi_), 0.5*M*xT*sin(phi_),
			      0.5*M*sqrt(max(sqr(x3)-sqr(xT),0.)),0.5*M*x3,ZERO);
    pphoton.rescaleEnergy();
    if(abs(pspect.z()+pemit.z()-pphoton.z())/M<1e-6) 
      pphoton.setZ(-pphoton.z());
    else if(abs(pspect.z()-pemit.z()+pphoton.z())/M<1e-6) 
      pemit .setZ(- pemit.z());
    // boost and rotate momenta
    LorentzRotation eventFrame( ( meMomenta()[2] + meMomenta()[3] ).findBoostToCM() );
    Lorentz5Momentum spectator;
    if(dipole==FF34) {
      spectator = eventFrame*meMomenta()[2];
    }
    else {
      spectator = eventFrame*meMomenta()[3];
    }
    eventFrame.rotateZ( -spectator.phi() );
    eventFrame.rotateY( -spectator.theta()  );
    eventFrame.invert();
    vector<Lorentz5Momentum> momenta(meMomenta());
    if(dipole==FF34) {
      momenta[3] = eventFrame*pspect;
      momenta[2] = eventFrame*pemit ;
    }
    else {
      momenta[2] = eventFrame*pspect;
      momenta[3] = eventFrame*pemit ;
    }
    momenta.push_back(eventFrame*pphoton);
    // phase-space factor
    double realFact = EMfact_*(1.-y)*jac*sqr(1.-2.*mu2)/sqrt(1.-4.*mu2);
    pair<double,double> realFF = make_pair(0.,0.);
    if(1.-x1>1e-5 && 1.-x2>1e-5) {
      realFF = subtractedQEDMEqqbar(momenta,dipole,true);
    }
    vector<double> output(2,0.);
    output[0] = realFF.first *realFact;
    output[1] = realFF.second*realFact;
    return output;
  }
  // ISR/FSR interference
  else if (dipole == IF13 || dipole == IF14 ||
	   dipole == IF23 || dipole == IF24) {
    // particles making up the dipole
    int iin   = dipole == IF13 || dipole == IF14 ? 0 : 1;
    int iout  = dipole == IF13 || dipole == IF23 ? 2 : 3;
    // momentum of system
    Lorentz5Momentum q = meMomenta()[iout]-meMomenta()[iin];
    q.rescaleMass();
    Energy  Q  = -q.mass();
    Energy2 Q2 = sqr(Q);
    // reduced mass
    Energy mj = meMomenta()[iout].mass();
    double mu2 = sqr(mj)/Q2;
    // second integration variable
    double vJac = (1.-mu2*z/(1.+mu2-z));
    double vt = vTilde_*vJac;
    // work out LT to Breit-Frame
    Lorentz5Momentum pb     = meMomenta()[iin ];
    Axis axis(q.vect().unit());
    double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
    LorentzRotation rot;
    if(axis.perp2()>1e-20) {
      rot.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
      rot.rotateX(Constants::pi);
    }
    if(abs(1.-q.e()/q.vect().mag())>1e-6) 
      rot.boostZ( q.e()/q.vect().mag());
    pb *= rot;
    if(pb.perp2()/GeV2>1e-20) {
      Boost trans = -1./pb.e()*pb.vect();
      trans.setZ(0.);
      rot.boost(trans);
    }
    // born incoming and outgoing momenta
    Lorentz5Momentum pij = rot*meMomenta()[iout];
    Lorentz5Momentum pat = rot*meMomenta()[iin ];
    // generate the real emission momenta
    double x1 = -(1.+mu2)/z;
    double x2 = 1-mu2-vt*(1.+mu2)/z;
    double x3 = 2.+x1-x2;
    double corr = (1.-z-vt)*mu2/((1.-vt)*(1.-z));
    double xT = sqrt(4.*(1.-z)/z*vt*(1.-vt)*(1.+corr));
    // incoming momentum
    Lorentz5Momentum pa(ZERO,ZERO,-0.5*x1*Q,-0.5*x1*Q,ZERO);
    // outgoing momentum
    Lorentz5Momentum pj( 0.5*Q*xT*cos(phi_), 0.5*Q*xT*sin(phi_),
			-0.5*Q*x2,0.5*Q*sqrt(sqr(x2)+sqr(xT)+4.*mu2),mj);
    Lorentz5Momentum pi(-0.5*Q*xT*cos(phi_),-0.5*Q*xT*sin(phi_),
			-0.5*Q*x3,0.5*Q*sqrt(sqr(x3)+sqr(xT)       ),ZERO);
    // boost back to the lab
    rot.invert();
    pa *= rot;
    pi *= rot;
    pj *= rot;
    // new momenta
    vector<Lorentz5Momentum> pnew(5);
    // incoming
    if(iin==0) {
      pnew[0] = pa;
      pnew[1] = meMomenta()[1];
    }
    else {
      pnew[0] = meMomenta()[0];
      pnew[1] = pa;
    }
    // outgoing
    if(iout==2) {
      pnew[2] = pj;
      pnew[3] = meMomenta()[3];
    }
    else {
      pnew[2] = meMomenta()[2];
      pnew[3] = pj;
    }
    // emitted
    pnew[4] = pi;
    // phase-space prefactors
    double phase = zJac*vJac/z*(1.+mu2);
    if      ( dipole == IF13 ) {
      realEmissionQEDIF13_ = pnew;
    }
    else if ( dipole == IF14 ) {
      realEmissionQEDIF14_ = pnew;
    }
    else if ( dipole == IF23 ) {
      realEmissionQEDIF23_ = pnew;
    }
    else if ( dipole == IF24 ) {
      realEmissionQEDIF24_ = pnew;
    }
    else
      assert(false);
    vector<double> output(2,0.);
    if(!(zTilde_<1e-7 || vt<1e-7 || 1.-z-vt < 1e-7 )) {
      pair<double,double> realIF = subtractedQEDMEqqbar(pnew,dipole,true);
      double fact = EMfact_*phase*newqPDF/oldqPDF;
      output[0] = realIF.first *fact;
      output[1] = realIF.second*fact;
    }
    // return the answer
    return output;
  }
  else {
    assert(false);
    return vector<double>(4,0.);
  }
}

pair<double,double> 
MEqq2gZ2ffPowhegQED::subtractedQEDMEqqbar(const vector<Lorentz5Momentum> & p,
					  DipoleType dipole, bool subtract) const {
  // calculate the matrix element
  cPDVector particles(mePartonData());
  particles.push_back(gamma_);
  vector<double> me = realQEDME(particles,p);
  // scale for the phase space
  Energy2 scale(ZERO);
  /////////// compute the two II dipole terms ////////////////////////////
  InvEnergy2 DII[2] = {ZERO,ZERO};
  if(QEDContributions_!=2) {
    //    cout<<"doing II dipoles"<<endl;
    double x = (p[0]*p[1]-p[4]*p[1]-p[4]*p[0])/(p[0]*p[1]);
    Lorentz5Momentum Kt = p[0]+p[1]-p[4];
    vector<Lorentz5Momentum> pa(4),pb(4);
    // momenta for q -> q gamma emission
    pa[0] = x*p[0];
    pa[1] =   p[1];
    Lorentz5Momentum K = pa[0]+pa[1];
    Lorentz5Momentum Ksum = K+Kt;
    Energy2 K2 = K.m2();
    Energy2 Ksum2 = Ksum.m2();
    for(unsigned int ix=2;ix<4;++ix)
      pa[ix] = p[ix]-2.*Ksum*(Ksum*p[ix])/Ksum2+2*K*(Kt*p[ix])/K2;
    // first LO matrix element
    double lo1 = loME(mePartonData(),pa,false);
    // momenta for qbar -> qbar gamma emission
    pb[0] =   p[0];
    pb[1] = x*p[1];
    K = pb[0]+pb[1];
    Ksum = K+Kt;
    K2 = K.m2();
    Ksum2 = Ksum.m2();
    for(unsigned int ix=2;ix<4;++ix)
      pb[ix] = p[ix]-2.*Ksum*(Ksum*p[ix])/Ksum2+2.*K*(Kt*p[ix])/K2;
    // second LO matrix element
    double lo2 = loME(mePartonData(),pb,false);
    // II dipoles
    DII[0] = 0.5/(p[0]*p[4])/x*(2./(1.-x)-(1.+x))*lo1;
    DII[1] = 0.5/(p[1]*p[4])/x*(2./(1.-x)-(1.+x))*lo2;
    // charge factors
    for(unsigned int ix=0;ix<2;++ix)
      DII[ix] *= sqr(double(mePartonData()[0]->iCharge())/3.);
    // phase-space factor
    if(dipole==II12||dipole==II21) scale = sHat();
  }
  /////////// compute the two FF dipole terms ////////////////////////////
  InvEnergy2 DFF[2]={ZERO,ZERO};
  Energy2 pT2final[2]={ZERO,ZERO};
  if(QEDContributions_!=1) {
    //cout<<"doing FF dipoles"<<endl;
    Lorentz5Momentum q = p[0]+p[1];
    Energy2 Q2=q.m2();
    Energy2 lambda = sqrt((Q2-sqr(p[2].mass()+p[3].mass()))*
			  (Q2-sqr(p[2].mass()-p[3].mass())));
    for(unsigned int iemit=0;iemit<2;++iemit) {
      unsigned int ispect = iemit==0 ? 1 : 0;
      Energy2 pipj = p[4      ] * p[2+iemit ];
      Energy2 pipk = p[4      ] * p[2+ispect];
      Energy2 pjpk = p[2+iemit] * p[2+ispect];
      double y = pipj/(pipj+pipk+pjpk);
      double z = pipk/(     pipk+pjpk);
      Energy mij = sqrt(2.*pipj+sqr(p[2+iemit].mass()));
      Energy2 lamB = sqrt((Q2-sqr(mij+p[2+ispect].mass()))*
			  (Q2-sqr(mij-p[2+ispect].mass())));
      Energy2 Qpk = q*p[2+ispect];
      Lorentz5Momentum pkt = 
	lambda/lamB*(p[2+ispect]-Qpk/Q2*q)
	+0.5/Q2*(Q2+sqr(p[2+ispect].mass())-sqr(p[2+ispect].mass()))*q;
      Lorentz5Momentum pijt = 
	q-pkt;
      double muj = p[2+iemit ].mass()/sqrt(Q2);
      double muk = p[2+ispect].mass()/sqrt(Q2);
      double vt = sqrt((1.-sqr(muj+muk))*(1.-sqr(muj-muk)))/(1.-sqr(muj)-sqr(muk));
      double v  = sqrt(sqr(2.*sqr(muk)+(1.-sqr(muj)-sqr(muk))*(1.-y))-4.*sqr(muk))
      /(1.-y)/(1.-sqr(muj)-sqr(muk));
      // transverse momentum
      double x1 = 2.*p[2+iemit ]*q/Q2;
      double x2 = 2.*p[2+ispect]*q/Q2;
      pT2final[iemit] = 0.25*Q2/(sqr(x2)-4.*sqr(muk))*
	((sqr(x1)-4.*sqr(muj))*(sqr(x2)-4.*sqr(muk)) - 
	 sqr(2.*x1+2.*x2-2.-x1*x2-2.*sqr(muj)-2.*sqr(muk)));
      // dipole term
      DFF[iemit] = 0.5/pipj*(2./(1.-(1.-z)*(1.-y))
			     -vt/v*(2.-z+sqr(p[2+iemit].mass())/pipj));
      // matrix element
      vector<Lorentz5Momentum> lomom(4);
      lomom[0] = p[0];
      lomom[1] = p[1];
      if(iemit==0) {
	lomom[2] = pijt;
	lomom[3] = pkt ;
      }
      else {
	lomom[3] = pijt;
	lomom[2] = pkt ;
      }
      // leading-order matrix element
      double lo = loME(mePartonData(),lomom,false);
      DFF[iemit] *= lo*sqr(double(mePartonData()[2]->iCharge())/3.);
    }
    // phase-space factpr
    if(dipole==FF34||dipole==FF43) scale = sHat();
  }
  /////////// ISR/FSR interference ////////////////////////////
  InvEnergy2 DIF[4]={ZERO,ZERO,ZERO,ZERO};
  Energy2 pT2IF[4] ={ZERO,ZERO,ZERO,ZERO};
  InvEnergy2 negativeDipoles(ZERO);
  if(QEDContributions_==0) {
//     cout<<"doing IF dipoles"<<endl;
    for(unsigned int ix=0;ix<4;++ix) {
      // incoming and outgoing particles in the dipole
      unsigned int iin  = ix<2 ? 0 :1;
      unsigned int iout = ix%2 == 0 ? 2 : 3;
      // emission variables
      double x = 1.-p[iout]*p[4]/(p[iin]*p[iout]+p[iin]*p[4]);
      double z =    p[iin ]*p[4]/(p[iin]*p[iout]+p[iin]*p[4]);
      // momenta for the born process
      Lorentz5Momentum pat = x*p[iin];
      Lorentz5Momentum pjt = p[4]+p[iout]-(1.-x)*p[iin];
      vector<Lorentz5Momentum> pa(4);
      for(unsigned iy=0;iy<4;++iy) {
	if(iy==iin)       pa[iy] = pat;
	else if(iy==iout) pa[iy] = pjt;
	else              pa[iy] = p[iy];
      }
      // pT
      Energy mj = p[iout].mass();
      Lorentz5Momentum q = pjt-pat;
      q.rescaleMass();
      Energy  Q  = -q.mass();
      Energy2 Q2 = sqr(Q);
      double mu2 = sqr(mj)/Q2;
      double corr = (1.-x-z)*mu2/((1.-z)*(1.-x));
      pT2IF[ix] = Q2*(1.-x)/x*z*(1.-z)*(1.+corr);
      // LO matrix element
      double lo = loME(mePartonData(),pa,false);
      Energy2 dot1 = max(p[iin ]*p[4],1e-30*MeV2);
      Energy2 dot2 = max(p[iout]*p[4],1e-30*MeV2);
      InvEnergy2 split = 
	0.5/dot1/x*(2./(1.-x+z)-1.-x) +
	0.5/dot2/x*(2./(1.-x+z)-2.+z-sqr(mj)/dot2);
      // lo piece and charge factors
      double charge = double(mePartonData()[iin]->iCharge()*
			     mePartonData()[iout]->iCharge())/9.;
      //new //ER
      // with the + sign because one of the 2 partons is outgoing,
      // so revert the sign of the
      // 'color' operator of the outgoing parton
      InvEnergy2 Dipole= + lo*charge*split; 
      if(DipoleSum_!=0 || charge>0)
	DIF[ix]          = Dipole;
      else if(dot1>1e-30*MeV2&&dot2>1e-30*MeV2)
	negativeDipoles += Dipole;
      if(int(ix)==int(dipole)-int(IF13)) scale = Q2;
    }
  }
  // supression function
  Energy2 pT2sup=ZERO;
  if(dipole==II12||dipole==II21) {
    pT2sup = sqr(p[4].x())+sqr(p[4].y());
  }
  else if(dipole==FF34||dipole==FF43) {
    pT2sup = dipole==FF34 ? pT2final[0] : pT2final[1];
  }
  else if(dipole==IF13) 
    pT2sup = pT2IF[0];
  else if(dipole==IF14) 
    pT2sup = pT2IF[1];
  else if(dipole==IF23) 
    pT2sup = pT2IF[2];
  else if(dipole==IF24) 
    pT2sup = pT2IF[3];
  else 
    assert(false);
  pair<double,double> supressionFactor = supressionFunction(pT2sup);
  // numerator for dipole ratio
  InvEnergy2 num;
  switch (dipole) {
  case II12 :
    num = DII[0];
    break;
  case II21 :
    num = DII[1];
    break;
  case FF34 :
    num = DFF[0];
    break;
  case FF43 :
    num = DFF[1];
    break;
  case IF13 :
    num = DIF[0];
    break;
  case IF14 :
    num = DIF[1];
    break;
  case IF23 :
    num = DIF[2];
    break;
  case IF24 :
    num = DIF[3];
    break;
  default:
    assert(false);
  };
  double meout(0.);
  InvEnergy2 den(ZERO);
  // ISR 
  if(QEDContributions_==1 ||
     (QEDContributions_==3 && (dipole==II12 || dipole ==II21))) {
    assert(dipole==II12 || dipole ==II21);
    meout = me[0];
    den = abs(DII[0])+abs(DII[1]);
    //    negativeDipoles=ZERO; //!ER:
  }
  // FSR 
  else if(QEDContributions_==2 ||
	  (QEDContributions_==3 && (dipole==FF34 || dipole ==FF43))) {
    assert(dipole==FF34 || dipole ==FF43);
    meout = me[1];
    den = abs(DFF[0])+abs(DFF[1]);
    //    negativeDipoles=ZERO; //!ER:
  }
  // full result
  else if(QEDContributions_==0) {
    meout = me[2];
    den = abs(DII[0])+abs(DII[1])+abs(DFF[0])+abs(DFF[1]);
    for(unsigned int ix=0;ix<4;++ix) den += abs(DIF[ix]);
  }
  else 
    assert(false);
  // return the answer
  pair<double,double> output = make_pair(0.,0.);
  if(den>ZERO) {
    if(subtract) {
      output.first  = scale*((UnitRemoval::InvE2*meout-negativeDipoles)*
			     abs(num)/den*supressionFactor.first -num);
    }
    else {
      output.first  = scale* UnitRemoval::InvE2*meout*
	abs(num)/den*supressionFactor.first;
    }
    output.second   = scale* UnitRemoval::InvE2*meout*
      abs(num)/den*supressionFactor.second;
  }
  return output;
}

pair<double,double>
MEqq2gZ2ffPowhegQED::subtractedQEDMEpqbar(const vector<Lorentz5Momentum> & p,
					  DipoleType dipole, bool subtract) const {
  // use the inheriting class to calculate the matrix element
  cPDVector particles(mePartonData());
  if(dipole==II12) {
     particles.push_back(particles[0]->CC());
     particles[0] = gluon_;
   }
  else if(dipole==II21) {
     particles.push_back(particles[1]->CC());
     particles[1] = gluon_;
   }
  else {
    assert(false);
  }
  vector<double> me = realQEDME(particles,p);
  // compute the two dipole terms
  double x = 1.-(p[4]*p[1]+p[4]*p[0])/(p[0]*p[1]);
  Lorentz5Momentum Kt = p[0]+p[1]-p[4];
  vector<Lorentz5Momentum> pa(4);
  // momenta for ISR
  if(dipole==II12) {
    pa[0] = x*p[0];
    pa[1] =   p[1];
  }
  else if(dipole==II21) {
    pa[0] =   p[0];
    pa[1] = x*p[1];
  }
  else {
    assert(false);
  }
  Lorentz5Momentum K = pa[0]+pa[1];
  Lorentz5Momentum Ksum = K+Kt;
  Energy2 K2 = K.m2();
  Energy2 Ksum2 = Ksum.m2();
  for(unsigned int ix=2;ix<4;++ix)
    pa[ix] = p[ix]-2.*Ksum*(Ksum*p[ix])/Ksum2+2*K*(Kt*p[ix])/K2; 
  // first LO matrix element 
  double lo1 = loME(mePartonData(),pa,false); 
  // dipole
  InvEnergy2 D1;
  if(dipole==II12) {
    D1 =  0.5/(p[0]*p[4])/x*(1.-2.*x*(1.-x));
  }
  else if(dipole==II21) {
    D1 =  0.5/(p[1]*p[4])/x*(1.-2.*x*(1.-x));
  }
  else {
    assert(false);
  }
  // charges
  D1 *= sqr(double(mePartonData()[2]->iCharge())/3.)*lo1;
  // supression function
  pair<double,double> supressionFactor = 
    supressionFunction(sqr(p[4].x())+sqr(p[4].y()));
  if(subtract) {
    return make_pair(sHat()*(UnitRemoval::InvE2*me[0]*supressionFactor.first -D1),
		     sHat()*(UnitRemoval::InvE2*me[0]*supressionFactor.second));
  }
  else {
    return make_pair(sHat()*(UnitRemoval::InvE2*me[0]*supressionFactor.first ),
		     sHat()*(UnitRemoval::InvE2*me[0]*supressionFactor.second));
  }
}

vector<double> 
MEqq2gZ2ffPowhegQED::realQEDME(const cPDVector & particles,
			       const vector<Lorentz5Momentum> & momenta) const {
  vector<SpinorWaveFunction>    fin(2),aout(2);
  vector<SpinorBarWaveFunction> ain(2),fout(2);
  vector<VectorWaveFunction> photon(2);
  // wavefunctions for the q qbar -> l+ l- gamma process
  if(particles[0]->id()==-particles[1]->id()) {
    for( unsigned int i = 0; i < 2; ++i ) {
      fin[i]    = SpinorWaveFunction   (momenta[0],particles[0],  i,incoming);
      ain[i]    = SpinorBarWaveFunction(momenta[1],particles[1],  i,incoming);
      photon[i] = VectorWaveFunction   (momenta[4],particles[4],2*i,outgoing);
    }
  }
  else if(particles[0]->id()==ParticleID::g && particles[1]->id()<0) {
    for( unsigned int i = 0; i < 2; ++i ) {
      fin[i]    = SpinorWaveFunction   (momenta[4],particles[4],  i,outgoing);
      ain[i]    = SpinorBarWaveFunction(momenta[1],particles[1],  i,incoming);
      photon[i] = VectorWaveFunction   (momenta[0],particles[0],2*i,incoming);
    }
  }
  else if(particles[0]->id()>0 && particles[1]->id()==ParticleID::g) {
    for( unsigned int i = 0; i < 2; ++i ) {
      fin[i]    = SpinorWaveFunction   (momenta[0],particles[0],  i,incoming);
      ain[i]    = SpinorBarWaveFunction(momenta[4],particles[4],  i,outgoing);
      photon[i] = VectorWaveFunction   (momenta[1],particles[1],2*i,incoming);
    }
  }
  else {
    for(unsigned int ix=0;ix<particles.size();++ix) {
      generator()->log() << particles[ix]->PDGName() << " " << momenta[ix]/GeV << "\n";
    }
    assert(false);
  }
  // wavefunctions for the leptons
  for(unsigned int i=0; i<2; ++i) {
    fout[i] = SpinorBarWaveFunction(momenta[2],particles[2],i,outgoing);
    aout[i] = SpinorWaveFunction   (momenta[3],particles[3],i,outgoing);
  }
  // off-shell vector bosons for speed
  Energy2 q2 = scale();
  VectorWaveFunction ZwaveIn[2][2],ZwaveOut[2][2];
  VectorWaveFunction PwaveIn[2][2],PwaveOut[2][2];
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      // intermediate Z
      ZwaveIn [ihel1][ihel2] = 
	oneLoopVertex_->evaluate(q2,1,Z0_   ,fin [ihel1],ain [ihel2]);
      ZwaveOut[ihel1][ihel2] = 
	oneLoopVertex_->evaluate(q2,1,Z0_   ,aout[ihel2],fout[ihel1]);
      // intermediate photon
      PwaveIn [ihel1][ihel2] = 
	oneLoopVertex_->evaluate(q2,1,gamma_,fin [ihel1],ain [ihel2]);
      PwaveOut[ihel1][ihel2] =
	oneLoopVertex_->evaluate(q2,1,gamma_,aout[ihel2],fout[ihel1]);
    }
  }
  vector<double> output(3,0.);
  vector<Complex> diagI(4,0.),diagF(4,0.);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int lhel1=0;lhel1<2;++lhel1) {
	for(unsigned int lhel2=0;lhel2<2;++lhel2) {
	  for(unsigned int ohel1=0;ohel1<2;++ohel1) {
	    // ISR diagrams
 	    // first Z diagram
 	    SpinorWaveFunction inters = oneLoopVertex_->
 	      evaluate(q2,5,fin[ihel1].particle(),fin[ihel1],photon[ohel1]);
  	    diagI[0] = oneLoopVertex_->
	      evaluate(q2, inters, ain[ihel2], ZwaveOut[lhel1][lhel2]);
 	    // second Z diagram
 	    SpinorBarWaveFunction interb = oneLoopVertex_->
 	      evaluate(q2,5,ain[ihel2].particle(),ain[ihel2],photon[ohel1]);
 	    diagI[1] = oneLoopVertex_->
	      evaluate(q2, fin[ihel1], interb, ZwaveOut[lhel1][lhel2]);
	    if(particles[2]->charged()) {
 	      // first photon diagram
	      diagI[2] = oneLoopVertex_->
		evaluate(q2, inters, ain[ihel2], PwaveOut[lhel1][lhel2]);
	      // second photon diagram
	      diagI[3] = oneLoopVertex_->
		evaluate(q2, fin[ihel1], interb, PwaveOut[lhel1][lhel2]);
	    }
	    // FSR diagrams
	    if(particles[2]->charged()) {
	      SpinorBarWaveFunction off1 = oneLoopVertex_->
		evaluate(q2,3,fout[lhel1].particle(),fout[lhel1],photon[ohel1]);
	      diagF[0] = oneLoopVertex_->
		evaluate(q2,aout[lhel2],off1,ZwaveIn[ihel1][ihel2]);
	      diagF[1] = oneLoopVertex_->
		evaluate(q2,aout[lhel2],off1,PwaveIn[ihel1][ihel2]);
	      SpinorWaveFunction    off2 = oneLoopVertex_->
		evaluate(q2,3,aout[lhel2].particle(),aout[lhel2],photon[ohel1]);
	      diagF[2] = oneLoopVertex_->
		evaluate(q2,off2,fout[lhel1],ZwaveIn[ihel1][ihel2]);
	      diagF[3] = oneLoopVertex_->
		evaluate(q2,off2,fout[lhel1],PwaveIn[ihel1][ihel2]);
	    }
	    // totals
	    Complex ISR = std::accumulate(diagI.begin(),diagI.end(),Complex(0.));
	    Complex FSR = std::accumulate(diagF.begin(),diagF.end(),Complex(0.));
	    output[0] += norm(ISR);
	    output[1] += norm(FSR);
	    output[2] += norm(ISR+FSR);
	  }
	}
      }
    }
  }
  // colour and spin factors
  if(particles[0]->id()==-particles[1]->id()) {
    for(unsigned int ix=0;ix<output.size();++ix)
      output[ix] *= 1./12.;
  }
  else  {
    for(unsigned int ix=0;ix<output.size();++ix)
      output[ix] *= 0.25;
  }
  // divided by 2 e^2
  for(unsigned int ix=0;ix<output.size();++ix)
    output[ix] *= 0.5/norm(oneLoopVertex_->norm());
  return output;
}

void MEqq2gZ2ffPowhegQED::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first );
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
  // order of particles
  if( hard[0]->id() != mePartonData()[0]->id() ) swap(hard[0], hard[1]);
  if( hard[2]->id() != mePartonData()[2]->id() ) swap(hard[2], hard[3]);
  // wavefunctions
  vector<SpinorWaveFunction>    fin,aout;
  vector<SpinorBarWaveFunction> ain,fout;
  SpinorWaveFunction(   fin ,hard[0],incoming,false,true);
  SpinorBarWaveFunction(ain ,hard[1],incoming,false,true);
  SpinorBarWaveFunction(fout,hard[2],outgoing,true ,true);
  SpinorWaveFunction(   aout,hard[3],outgoing,true ,true);
  qqbarME(fin,ain,fout,aout,true);
  // get the spin info objects
  SpinPtr spin[4];
  for(unsigned int ix=0;ix<4;++ix) spin[ix]=hard[ix]->spinInfo();
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(me_);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix) 
    spin[ix]->productionVertex(hardvertex);
}

void MEqq2gZ2ffPowhegQED::hardQCDEmission(vector<ShowerProgenitorPtr> & 
					  particlesToShower, int & emission_type,
					  Energy & pTmax) {
  emission_type = -1;
  pTmax = -GeV;
  Energy rootS = sqrt(lastS());
  // limits on the rapidity of the jet
  double minyj = -10.0,maxyj = 10.0;
  pair<double,double> x = make_pair(particlesToShower[0]->progenitor()->x(),
				    particlesToShower[1]->progenitor()->x());
  // loop over the possible emissions
  unsigned int imin(0),imax(4);
  if(QCDRadiationType_>=0) {
    imin = QCDRadiationType_;
    imax = imin+1;
  }
  for(unsigned int ix=imin;ix<imax;++ix) {
    Energy pT = 0.5*generator()->maximumCMEnergy();
    // particles for the hard process
    cPDVector particles;
    for(unsigned int iy=0;iy<particlesToShower.size();++iy) {
      particles.push_back(particlesToShower[iy]->progenitor()->dataPtr());
    }
    if(ix<2) particles.push_back(gluon_);
    else if(ix==2) {
      particles.push_back(particles[0]->CC());
      particles[0] = gluon_;
    }
    else {
      particles.push_back(particles[1]->CC());
      particles[1] = gluon_;
    }
    vector<Lorentz5Momentum> momenta(5);
    double a = alphaQCD_->overestimateValue()/Constants::twopi*
      prefactorQCD_[ix]*(maxyj-minyj);
    do {
      pT *= pow(UseRandom::rnd(),1./a);
      if(pT<=minpTQCD_) break;
      double y = UseRandom::rnd()*(maxyj-minyj)+ minyj;
      double vt,z;
      if(ix%2==0) {
	vt = pT*exp(-y)/rootS/x.second;
	z  = (1.-pT*exp(-y)/rootS/x.second)/(1.+pT*exp( y)/rootS/x.first );
	if(z>1.||z<x.first) continue;
      }
      else {
	vt = pT*exp( y)/rootS/x.first ;
	z  = (1.-pT*exp( y)/rootS/x.first )/(1.+pT*exp(-y)/rootS/x.second );
	if(z>1.||z<x.second) continue;
      }
      if(vt>1.-z || vt<0.) continue;
      if(ix%2==0) {
	momenta[0] = particlesToShower[0]->progenitor()->momentum()/z;
	momenta[1] = particlesToShower[1]->progenitor()->momentum();
      }
      else {
	momenta[0] = particlesToShower[0]->progenitor()->momentum();
	momenta[1] = particlesToShower[1]->progenitor()->momentum()/z;
      }
      double phi = Constants::twopi*UseRandom::rnd();
      momenta[2] = particlesToShower[2]->progenitor()->momentum();
      momenta[3] = particlesToShower[3]->progenitor()->momentum();
      if(!_quarkplus) y *= -1.;
      momenta[4] = Lorentz5Momentum(pT*cos(phi),pT*sin(phi),
				    pT*sinh(y),pT*cosh(y), ZERO);
      Lorentz5Momentum K = momenta[0] + momenta[1] - momenta[4]; 
      Lorentz5Momentum Kt = momenta[2]+momenta[3];
      Lorentz5Momentum Ksum = K+Kt;
      Energy2 K2 = Kt.m2(), Ksum2 = Ksum.m2();
      for(unsigned int iy=2;iy<4;++iy) {
	momenta [iy] = momenta [iy] - 2.*Ksum*(Ksum*momenta [iy])/Ksum2
	  +2.*K*(Kt*momenta [iy])/K2;
      }
      // matrix element piece
      double wgt = alphaQCD_->ratio(sqr(pT))*z/(1.-vt)/prefactorQCD_[ix]/loME_;
      // compute me piece here
      if(ix==0)
	wgt *= 8./3.*sqr(pT)/sHat()*subtractedQCDMEqqbar(momenta,II12,false).first;
      else if(ix==1)
	wgt *= 8./3.*sqr(pT)/sHat()*subtractedQCDMEqqbar(momenta,II21,false).first;
      else if(ix==2)
	wgt *=       sqr(pT)/sHat()*subtractedQCDMEgqbar(momenta,II12,false).first;
      else if(ix==3)
 	wgt *=       sqr(pT)/sHat()*subtractedQCDMEgqbar(momenta,II21,false).first;
      // pdf piece
      double pdf[4];
      if(ix%2==0) {
	pdf[0] = _beams[0]->pdf()->xfx(_beams[0],_partons [0],
				       scale()        ,x.first   )  /x.first;
	pdf[1] = _beams[0]->pdf()->xfx(_beams[0],particles[0],
				       scale()+sqr(pT),x.first /z)*z/x.first;
	pdf[2] = _beams[1]->pdf()->xfx(_beams[1],_partons [1],
				       scale()        ,x.second  )  /x.second;
	pdf[3] = _beams[1]->pdf()->xfx(_beams[1],particles[1],
				       scale()+sqr(pT),x.second  )  /x.second;
      }
      else {
	pdf[0] = _beams[1]->pdf()->xfx(_beams[1],_partons [1],
				       scale()        ,x.second  )  /x.second;
	pdf[1] = _beams[1]->pdf()->xfx(_beams[1],particles[1],
				       scale()+sqr(pT),x.second/z)*z/x.second;
	pdf[2] = _beams[0]->pdf()->xfx(_beams[0],_partons [0],
				       scale()        ,x.first   )  /x.first;
	pdf[3] = _beams[0]->pdf()->xfx(_beams[0],particles[0],
				       scale()+sqr(pT),x.first   )  /x.first;
      }
      if(pdf[0]<=0.||pdf[1]<=0.) continue;
      if(pdf[2]<=0.||pdf[3]<=0.) continue;
      wgt *= pdf[1]/pdf[0];
      wgt *= pdf[3]/pdf[2];
      // check weight less than one
      if(wgt>1.) {
	generator()->log() << "Weight greater than one for emission type " << ix
			   << "in MEqq2gZ2ffPowhegQED::generateHardest()"
			   << " weight = " << wgt << "\n";
      }
      // break if select emission
      if(UseRandom::rnd()<wgt) break;
    }
    while(pT>minpTQCD_);
    if(pT>minpTQCD_ && pT>pTmax) {
      pTmax = pT;
      if(ix==0) {
	realEmissionQCDGluon1_=momenta;
	emission_type=1;
      }
      else if(ix==1) {
	realEmissionQCDGluon2_=momenta;
	emission_type=3;
      }
      else if(ix==2) {
	realEmissionQCDQuark1_=momenta;
	emission_type=2;
      }
      else if(ix==3) {
	realEmissionQCDQuark2_=momenta;
	emission_type=4;
      }
    }
  }
}

/**
 *  Generate a hard QED emission
 */
void MEqq2gZ2ffPowhegQED::hardQEDEmission(vector<ShowerProgenitorPtr> & particlesToShower,
					 int & emission_type, Energy & pTmax) {
  emission_type = -1;
  pTmax = -GeV;
  Energy pT_II(-GeV),pT_IF(-GeV),pT_FF(-GeV);
  int type_II(-1),type_IF(-1),type_FF(-1);
  // initial-initial radiation
  hardQEDIIEmission(particlesToShower,type_II,pT_II);
  // final-final radiation
  hardQEDFFEmission(particlesToShower,type_FF,pT_FF);
  // initial-final radiation
  hardQEDIFEmission(particlesToShower,type_IF,pT_IF);
  if(type_II<0&&type_FF<0&&type_IF<0) return;
  // select the maximum pT emission
  if( pT_II>pT_FF && pT_II>pT_IF) {
    pTmax = pT_II;
    emission_type = type_II;
  }
  else if (pT_FF>pT_II && pT_FF>pT_IF) {
    pTmax = pT_FF;
    emission_type = type_FF;
  }
  else if (pT_IF>pT_FF && pT_IF>pT_II) {
    pTmax = pT_IF;
    emission_type = type_IF;
  }
}

void MEqq2gZ2ffPowhegQED::hardQEDIIEmission(vector<ShowerProgenitorPtr> & particlesToShower,
					    int & emission_type, Energy & pTmax) {
  emission_type = -1;
  pTmax = -GeV;
  if(QEDContributions_==2) return;
  Energy rootS = sqrt(lastS());
  // limits on the rapidity of the jet
  double minyj = -10.0,maxyj = 10.0;
  pair<double,double> x = make_pair(particlesToShower[0]->progenitor()->x(),
				    particlesToShower[1]->progenitor()->x());
  // loop over the possible emissions
  vector<Energy> pT;
  double charge = sqr(double(particlesToShower[0]->progenitor()->dataPtr()->iCharge())/3.);
  for(unsigned int ix=0;ix<4;++ix) {
    // skip incoming photons if not required
    if(!incomingPhotons_&&(ix==2||ix==3)) continue;
    pT.push_back(0.5*generator()->maximumCMEnergy());
    // particles for the hard process
    cPDVector particles;
    for(unsigned int iy=0;iy<particlesToShower.size();++iy) {
      particles.push_back(particlesToShower[iy]->progenitor()->dataPtr());
    }
    if(ix<2) particles.push_back(gamma_);
    else if(ix==2) {
      particles.push_back(particles[0]->CC());
      particles[0] = gamma_;
    }
    else {
      particles.push_back(particles[1]->CC());
      particles[1] = gamma_;
    }
    vector<Lorentz5Momentum> momenta(5);
    double a = alphaQED_->overestimateValue()/Constants::twopi*
      prefactorQED_[ix]*(maxyj-minyj)*charge;
    do {
      pT[ix] *= pow(UseRandom::rnd(),1./a);
      if(pT[ix]<=minpTQEDII_) break;
      double y = UseRandom::rnd()*(maxyj-minyj)+ minyj;
      double vt,z;
      if(ix%2==0) {
	vt = pT[ix]*exp(-y)/rootS/x.second;
	z  = (1.-pT[ix]*exp(-y)/rootS/x.second)/(1.+pT[ix]*exp( y)/rootS/x.first );
	if(z>1.||z<x.first) continue;
      }
      else {
	vt = pT[ix]*exp( y)/rootS/x.first ;
	z  = (1.-pT[ix]*exp( y)/rootS/x.first )/(1.+pT[ix]*exp(-y)/rootS/x.second );
	if(z>1.||z<x.second) continue;
      }
      if(vt>1.-z || vt<0.) continue;
      if(ix%2==0) {
	momenta[0] = particlesToShower[0]->progenitor()->momentum()/z;
	momenta[1] = particlesToShower[1]->progenitor()->momentum();
      }
      else {
	momenta[0] = particlesToShower[0]->progenitor()->momentum();
	momenta[1] = particlesToShower[1]->progenitor()->momentum()/z;
      }
      double phi = Constants::twopi*UseRandom::rnd();
      momenta[2] = particlesToShower[2]->progenitor()->momentum();
      momenta[3] = particlesToShower[3]->progenitor()->momentum();
      if(!_quarkplus) y *= -1.;
      momenta[4] = Lorentz5Momentum(pT[ix]*cos(phi),pT[ix]*sin(phi),
				    pT[ix]*sinh(y),pT[ix]*cosh(y), ZERO);
      Lorentz5Momentum K = momenta[0] + momenta[1] - momenta[4]; 
      Lorentz5Momentum Kt = momenta[2]+momenta[3];
      Lorentz5Momentum Ksum = K+Kt;
      Energy2 K2 = K.m2(), Ksum2 = Ksum.m2();
      for(unsigned int iy=2;iy<4;++iy) {
	momenta [iy] = momenta [iy] - 2.*Ksum*(Ksum*momenta [iy])/Ksum2
	  +2*K*(Kt*momenta [iy])/K2;
      }
      // matrix element piece
      double wgt = alphaQED_->ratio(sqr(pT[ix]))*z/(1.-vt)/prefactorQED_[ix]/loME_/charge;
      // compute me piece here
      if(ix==0)
	wgt *= 2.*sqr(pT[ix])/sHat()*subtractedQEDMEqqbar(momenta,II12,false).first;
      else if(ix==1)
	wgt *= 2.*sqr(pT[ix])/sHat()*subtractedQEDMEqqbar(momenta,II21,false).first;
      else if(ix==2)
	wgt *= 2.*sqr(pT[ix])/sHat()*subtractedQEDMEpqbar(momenta,II12,false).first;
      else if(ix==3)
	wgt *= 2.*sqr(pT[ix])/sHat()*subtractedQEDMEpqbar(momenta,II21,false).first;
      // pdf piece
      double pdf[2];
      if(ix%2==0) {
	pdf[0] = _beams[0]->pdf()->xfx(_beams[0],_partons [0],
				       scale(),            x.first   )  /x.first;
	pdf[1] = _beams[0]->pdf()->xfx(_beams[0],particles[0],
				       scale()+sqr(pT[ix]),x.first /z)*z/x.first;
      }
      else {
	pdf[0] = _beams[1]->pdf()->xfx(_beams[1],_partons [1],
				       scale()            ,x.second  )  /x.second;
	pdf[1] = _beams[1]->pdf()->xfx(_beams[1],particles[1],
				       scale()+sqr(pT[ix]),x.second/z)*z/x.second;
      }
      if(pdf[0]<=0.||pdf[1]<=0.) continue;
      wgt *= pdf[1]/pdf[0];
      // check weight less than one
      if(wgt>1.) {
	generator()->log() << "Weight greater than one for emission type " << ix
			   << "in MEqq2gZ2ffPowhegQED::hardQEDIIEmission()"
			   << " weight = " << wgt << "\n";
      }
      // break if select emission
      if(UseRandom::rnd()<wgt) break;
    }
    while(pT[ix]>minpTQEDII_);
    if(pT[ix]>minpTQEDII_ && pT[ix]>pTmax) {
      pTmax = pT[ix];
      if(ix==0) {
	realEmissionQEDPhoton1_ = momenta;
	emission_type=5;
      }
      else if(ix==1) {
	realEmissionQEDPhoton2_ = momenta;
	emission_type=7;
      }
      else if(ix==2) {
	realEmissionQEDQuark1_  = momenta;
	emission_type=6;
      }
      else if(ix==3) {
	realEmissionQEDQuark2_  = momenta;
	emission_type=8;
      }
    }
  }
}

void MEqq2gZ2ffPowhegQED::hardQEDIFEmission(vector<ShowerProgenitorPtr> & particlesToShower,
					    int & emission_type, Energy & pTmax) {
  if(QEDContributions_!=0) return;
  emission_type = -1;
  pTmax = -GeV;
  // extract the particles
  cPDVector particles;
  for(unsigned int iy=0;iy<particlesToShower.size();++iy) {
    particles.push_back(particlesToShower[iy]->progenitor()->dataPtr());
  }
  particles.push_back(gamma_);
  // momentum fractions of the incoming particles
  pair<double,double> x = make_pair(particlesToShower[0]->progenitor()->x(),
				    particlesToShower[1]->progenitor()->x());
  // loop over the possible dipoles
  for(unsigned int idipole=0;idipole<4;++idipole) {
    // incoming and outgoing particles in the dipole
    unsigned int iin  = idipole<2 ? 0 :1;
    unsigned int iout = idipole%2 == 0 ? 2 : 3;
    // charge of the dipole
    double charge = +
      particlesToShower[iin ]->progenitor()->dataPtr()->iCharge()*
      particlesToShower[iout]->progenitor()->dataPtr()->iCharge()/9.;
    if(charge<0.) continue;
    // which dipole are we doing
    DipoleType dipole;
    if     (idipole==0) dipole = IF13;
    else if(idipole==1) dipole = IF14;
    else if(idipole==2) dipole = IF23;
    else if(idipole==3) dipole = IF24;
    // storage of the real emission momenta
    vector<Lorentz5Momentum> realMomenta(5);
    Lorentz5Momentum q = 
      particlesToShower[iout]->progenitor()->momentum() - 
      particlesToShower[iin ]->progenitor()->momentum();
    q.rescaleMass();
    Energy  Q  = -q.mass();
    Energy2 Q2 = sqr(Q);
    double xB = iin==0 ? x.first : x.second;
    // reduced mass
    Energy mj = particlesToShower[iout]->progenitor()->momentum().mass();
    double mu2 = sqr(mj)/Q2;
    // work out LT to the Breit-Frame
    Lorentz5Momentum pb     = particlesToShower[iin ]->progenitor()->momentum();
    Axis axis(q.vect().unit());
    double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
    LorentzRotation rot;
    if(axis.perp2()>1e-20) {
      rot.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
      rot.rotateX(Constants::pi);
    }
    if(abs(1.-q.e()/q.vect().mag())>1e-6) 
      rot.boostZ( q.e()/q.vect().mag());
    pb *= rot;
    if(pb.perp2()/GeV2>1e-20) {
      Boost trans = -1./pb.e()*pb.vect();
      trans.setZ(0.);
      rot.boost(trans);
    }
    rot.invert();
    // momenta not effected by the emission
    for(unsigned int iy=0;iy<4;++iy) {
      if(iy==iin||iy==iout) continue;
      realMomenta[iy] = particlesToShower[iy]->progenitor()->momentum();
    }
    // maximum value of the xT
    double xT = sqrt((1.-xB)/xB);
    double xTMin = 2.*minpTQEDIF_/sqrt(Q2);
    double wgt(0.);
    // prefactor
    double a = charge*alphaQED_->overestimateValue()*preIFQED_/Constants::twopi;
    // loop to generate kinematics
    double xp,zp;
    double n=0.;
    do {
      wgt = 0.;
      // intergration variables dxT/xT^3
      xT *= pow(1.-(n+2.)*log(UseRandom::rnd())/a*pow(xT,n+2.),-1./(n+2));
//       xT *= 1./sqrt(1.-2.*log(UseRandom::rnd())/a*sqr(xT));
      // zp
      zp = UseRandom::rnd();
      // xp
      xp = (1-zp*mu2/(1.-zp+mu2))/(1.+0.25*sqr(xT)/zp/(1.-zp+mu2));
      // check allowed
      if(xp<xB||xp>1.) continue;
      // azimuth
      double phi = Constants::twopi*UseRandom::rnd();
      double x1 = -(1.+mu2)/xp;
      double x2 = 1.-mu2-zp*(1.+mu2)/xp;
      double x3 = 2.+x1-x2;
      // now compute the momenta
      realMomenta[iin] = 
	Lorentz5Momentum(ZERO,ZERO,-0.5*x1*Q,-0.5*x1*Q,ZERO);
      realMomenta[iout] = 
	Lorentz5Momentum( 0.5*Q*xT*cos(phi), 0.5*Q*xT*sin(phi),
			 -0.5*Q*x2,0.5*Q*sqrt(sqr(x2)+sqr(xT)+4.*mu2),mj);
      realMomenta[4] = 
	Lorentz5Momentum(-0.5*Q*xT*cos(phi),-0.5*Q*xT*sin(phi),
			 -0.5*Q*x3,0.5*Q*sqrt(sqr(x3)+sqr(xT)),ZERO);
      realMomenta[iin]  *= rot;
      realMomenta[iout] *= rot;
      realMomenta[4]    *= rot;
      // phase-space piece of the weight
      double corr =(1.-xp-zp)*mu2/(1.-xp)/(1.-zp);
      wgt = pow(xT,n)*8.*zp*(1.-zp)*sqr(1.-xp)/xp/preIFQED_*(1.+mu2)/xp*
	sqr(1.+corr)/(1.+mu2);
      wgt *= alphaQED_->ratio(0.25*Q2*sqr(xT));
      // PDF piece
      double pdf[2]={_beams[iin]->pdf()->xfx(_beams[iin],_partons [iin],
					     scale(),xB)  /xB,
		     _beams[iin]->pdf()->xfx(_beams[iin],particles[iin],
					     0.5*xT*Q2,xB/xp)*xp/xB};
      
      if(pdf[0]<=0.||pdf[1]<=0.) {
	wgt =0.;
	continue;
      }
      wgt *= pdf[1]/pdf[0];
      // matrix element piece
      wgt *= subtractedQEDMEqqbar(realMomenta,dipole,false).first/charge/loME_;
      if(wgt>1.) {
	generator()->log() << "problem with IF weight " << wgt << " " << idipole << " "
			   << xT << " " << zp << " " << xp << "\n";
	for(unsigned int ix=0;ix<realMomenta.size();++ix)
	  generator()->log() << ix << " " << realMomenta[ix]/GeV << "\n";
	generator()->log() << "testing masses "
			   << (realMomenta[0]+realMomenta[1]).m()/GeV << " "
			   << (realMomenta[2]+realMomenta[3]).m()/GeV << "\n";
	generator()->log() << "sum "
			   << (realMomenta[0]+realMomenta[1]-realMomenta[2]-realMomenta[3]-realMomenta[4])/GeV << "\n";
      }
    }
    while(xT>xTMin&&UseRandom::rnd()>wgt);
    if(xT<=xTMin) continue;
    Energy pT = xT*0.5*Q;
    if(pT<pTmax) continue;
    pTmax = pT;
    int isif = zp>xp ? 4 : 0;
    emission_type = 11+idipole+isif;
    if(idipole==0) 
      realEmissionQEDIF13_ = realMomenta;
    else if(idipole==1)
      realEmissionQEDIF14_ = realMomenta;
    else if(idipole==2)
      realEmissionQEDIF23_ = realMomenta;
    else if(idipole==3)
      realEmissionQEDIF24_ = realMomenta;
  }
}

void MEqq2gZ2ffPowhegQED::
hardQEDFFEmission(vector<ShowerProgenitorPtr> & particlesToShower,
		  int & emission_type, Energy & pTmax) {
  if(QEDContributions_==1) return;
  emission_type = -1;
  pTmax = -GeV;
  // extract the particles
  cPDVector particles;
  for(unsigned int iy=0;iy<particlesToShower.size();++iy) {
    particles.push_back(particlesToShower[iy]->progenitor()->dataPtr());
  }
  particles.push_back(gamma_);
  // boost from lab to CMS frame with outgoing particles
  // along the z axis
  Lorentz5Momentum pcms = particlesToShower[2]->progenitor()->momentum() + 
                          particlesToShower[3]->progenitor()->momentum();
  LorentzRotation eventFrame( pcms.findBoostToCM() );
  Lorentz5Momentum spectator = eventFrame*particlesToShower[2]->progenitor()->momentum();
  eventFrame.rotateZ( -spectator.phi() );
  eventFrame.rotateY( -spectator.theta()  );
  eventFrame.invert();
  // mass of the final-state system
  Energy2 M2 = pcms.m2();
  Energy  M  = sqrt(M2);
  double mu1 = particlesToShower[2]->progenitor()->momentum().mass()/M;
  double mu2 = particlesToShower[2]->progenitor()->momentum().mass()/M;
  double mu12 = sqr(mu1), mu22 = sqr(mu2);
  double lambda = sqrt(1.+sqr(mu12)+sqr(mu22)-2.*mu12-2.*mu22-2.*mu12*mu22);
  // max pT
  pTmax = 0.5*sqrt(M2)*(1.-sqr(mu1+mu2));
  // max y
  double ymax = acosh(pTmax/minpTQEDFF_);
  if(!particlesToShower[2]->progenitor()->dataPtr()->charged())
    return;
  double charge = sqr(double(particlesToShower[2]->
			     progenitor()->dataPtr()->iCharge())/3.);
  double a = alphaQED_->overestimateValue()/Constants::twopi*2.*ymax*preFFQED_*charge;
  // variables for the emission
  Energy pT[2];
  double  y[2],phi[2],x3[2],x1[2][2],x2[2][2];
  double contrib[2][2];
  // storage of the real emission momenta
  vector<Lorentz5Momentum> realMomenta[2][2]=
    {{vector<Lorentz5Momentum>(5),vector<Lorentz5Momentum>(5)},
     {vector<Lorentz5Momentum>(5),vector<Lorentz5Momentum>(5)}};
  for(unsigned int ix=0;ix<2;++ix)
    for(unsigned int iy=0;iy<2;++iy)
      for(unsigned int iz=0;iz<2;++iz)
	realMomenta[ix][iy][iz] = particlesToShower[iz]->progenitor()->momentum();
  // generate the emission
  for(unsigned int ix=0;ix<2;++ix) {
    if(ix==1) {
      swap(mu1 ,mu2 );
      swap(mu12,mu22);
    }
    pT[ix] = pTmax;
    y [ix] = 0.;
    bool reject = true;
    do {
      // generate pT
      pT[ix] *= pow(UseRandom::rnd(),1./a);
      if(pT[ix]<minpTQEDFF_) {
	pT[ix] = -GeV;
	break;
      }
      // generate y
      y[ix] = -ymax+2.*UseRandom::rnd()*ymax;
      // generate phi
      phi[ix] = UseRandom::rnd()*Constants::twopi;
      // calculate x3 and check in allowed region
      x3[ix] = 2.*pT[ix]*cosh(y[ix])/M;
      if(x3[ix] < 0. || x3[ix] > 1. -sqr( mu1 + mu2 ) ) continue;
      // find the possible solutions for x1
      double xT2 = sqr(2./M*pT[ix]);
      double root = (-sqr(x3[ix])+xT2)*
	(xT2*mu22+2.*x3[ix]-sqr(mu12)+2.*mu22+2.*mu12-sqr(x3[ix])-1.
	 +2.*mu12*mu22-sqr(mu22)-2.*mu22*x3[ix]-2.*mu12*x3[ix]);
      double c1=2.*sqr(x3[ix])-4.*mu22-6.*x3[ix]+4.*mu12-xT2*x3[ix]
	+2.*xT2-2.*mu12*x3[ix]+2.*mu22*x3[ix]+4.;
      if(root<0.) continue;
      x1[ix][0] = 1./(4.-4.*x3[ix]+xT2)*(c1-2.*sqrt(root));
      x1[ix][1] = 1./(4.-4.*x3[ix]+xT2)*(c1+2.*sqrt(root));
      // change sign of y if 2nd particle emits
      if(ix==1) y[ix] *=-1.;
      // loop over the solutions
      for(unsigned int iy=0;iy<2;++iy) {
	contrib[ix][iy]=0.;
	// check x1 value allowed
	if(x1[ix][iy]<2.*mu1||x1[ix][iy]>1.+mu12-mu22) continue;
	// calculate x2 value and check allowed
	x2[ix][iy] = 2.-x3[ix]-x1[ix][iy];
	double root = max(0.,sqr(x1[ix][iy])-4.*mu12);
	root = sqrt(root);
	double x2min = 1.+mu22-mu12
	  -0.5*(1.-x1[ix][iy]+mu12-mu22)/(1.-x1[ix][iy]+mu12)*(x1[ix][iy]-2.*mu12+root);
	double x2max = 1.+mu22-mu12
	  -0.5*(1.-x1[ix][iy]+mu12-mu22)/(1.-x1[ix][iy]+mu12)*(x1[ix][iy]-2.*mu12-root);
	if(x2[ix][iy]<x2min||x2[ix][iy]>x2max) continue;
	// check the z components
	double z1 =  sqrt(sqr(x1[ix][iy])-4.*mu12-xT2);
	double z2 = -sqrt(sqr(x2[ix][iy])-4.*mu22);
	double z3 =  pT[ix]*sinh(y[ix])*2./M;
	if(ix==1) z3 *=-1.;
	if(abs(-z1+z2+z3)<1e-9) z1 *= -1.;
	if(abs(z1+z2+z3)>1e-5) continue;
	// construct the momenta
	realMomenta[ix][iy][4] =
	  Lorentz5Momentum(pT[ix]*cos(phi[ix]),pT[ix]*sin(phi[ix]),
			   pT[ix]*sinh(y[ix]) ,pT[ix]*cosh(y[ix]),ZERO);
	if(ix==0) {
	  realMomenta[ix][iy][2] =
	    Lorentz5Momentum(-pT[ix]*cos(phi[ix]),-pT[ix]*sin(phi[ix]),
			     z1*0.5*M,x1[ix][iy]*0.5*M,M*mu1);
	  realMomenta[ix][iy][3] =
	    Lorentz5Momentum(ZERO,ZERO, z2*0.5*M,x2[ix][iy]*0.5*M,M*mu2);
	}
	else {
	  realMomenta[ix][iy][2] =
	    Lorentz5Momentum(ZERO,ZERO,-z2*0.5*M,x2[ix][iy]*0.5*M,M*mu2);
	  realMomenta[ix][iy][3] =
	    Lorentz5Momentum(-pT[ix]*cos(phi[ix]),-pT[ix]*sin(phi[ix]),
			     -z1*0.5*M,x1[ix][iy]*0.5*M,M*mu1);
	}
	// boost the momenta back to the lab
	for(unsigned int iz=2;iz<5;++iz)
	  realMomenta[ix][iy][iz] *= eventFrame;
	// jacobian and prefactors for the weight
	Energy J = M/sqrt(xT2)*abs(-x1[ix][iy]*x2[ix][iy]+2.*mu22*x1[ix][iy]
				   +x2[ix][iy]+x2[ix][iy]*mu12+mu22*x2[ix][iy]
				   -sqr(x2[ix][iy]))
	  /pow(sqr(x2[ix][iy])-4.*mu22,1.5);
	// prefactors etc
	contrib[ix][iy] = 0.5*pT[ix]/J/preFFQED_/lambda;
	// matrix element piece
	double me;
	if(ix==0) me = subtractedQEDMEqqbar(realMomenta[ix][iy],
					    FF34,false).first/charge/loME_;
	else      me = subtractedQEDMEqqbar(realMomenta[ix][iy],
					    FF43,false).first/charge/loME_;
	contrib[ix][iy] *= 2.*me;
	// coupling piece
	contrib[ix][iy] *= alphaQED_->ratio(sqr(pT[ix]));
      }
      if(contrib[ix][0]+contrib[ix][1]>1.) {
	ostringstream s;
	s << "MEee2gZ2qq::generateHardest weight for channel " << ix
	  << "is " << contrib[ix][0]+contrib[ix][1] 
	  << " which is greater than 1";
	generator()->logWarning( Exception(s.str(), Exception::warning) );
      }
      reject =  UseRandom::rnd() > contrib[ix][0] + contrib[ix][1];
    }
    while (reject);
    if(pT[ix]<minpTQEDFF_)
      pT[ix] = -GeV;
  }
  emission_type = -1;
  pTmax = -GeV;
  if(pT[0]<ZERO&&pT[1]<ZERO) return;
  // now pick the emmision with highest pT
  if(pT[0]>pT[1]) {
    emission_type=9;
    pTmax = pT[0];
    if(UseRandom::rnd()<contrib[0][0]/(contrib[0][0]+contrib[0][1]))
      realEmissionQEDFFLepton3_ = realMomenta[0][0];
    else
      realEmissionQEDFFLepton3_ = realMomenta[0][1];
  }
  else {
    emission_type=10;
    pTmax = pT[1];
    if(UseRandom::rnd()<contrib[1][0]/(contrib[1][0]+contrib[1][1]))
      realEmissionQEDFFLepton4_ = realMomenta[1][0];
    else
      realEmissionQEDFFLepton4_ = realMomenta[1][1];
  }
}

HardTreePtr MEqq2gZ2ffPowhegQED::generateHardest(ShowerTreePtr tree,
						 vector<ShowerInteraction::Type> inter) {
  // get the particles to be showered
  _beams.clear();
  _partons.clear();
  // find the incoming particles
  ShowerParticleVector incoming;
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  _quarkplus = true;
  vector<ShowerProgenitorPtr> particlesToShower;
  //progenitor particles are produced in z direction.
  for( cit = tree->incomingLines().begin(); cit != tree->incomingLines().end(); ++cit ) {
    incoming.push_back( cit->first->progenitor() );
    _beams.push_back( cit->first->beam() );
    _partons.push_back( cit->first->progenitor()->dataPtr() );
    // check that quark is along +ve z direction
    if(cit->first->progenitor()->id() > 0 &&
       cit->first->progenitor()->momentum().z() < ZERO ) 
      _quarkplus = false;
    particlesToShower.push_back( cit->first );
  }
  // we are assuming quark first, swap order to ensure this
  // if antiquark first
  if(_partons[0]->id()<_partons[1]->id()) {
    swap(_partons[0],_partons[1]);
    swap(_beams[0],_beams[1]);
    swap(incoming[0],incoming[1]);
    swap(particlesToShower[0],particlesToShower[1]);
  }
  // outgoing particles
  for( map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	 cjt= tree->outgoingLines().begin();
       cjt != tree->outgoingLines().end();++cjt ) {
    particlesToShower.push_back( cjt->first );
  }
  if(particlesToShower[2]->id()<0)
    swap(particlesToShower[2],particlesToShower[3]);
  // genuine hard emission in matrix element
  double wgtb = std::accumulate(++weights_.begin(),weights_.end(),0.);
  int emission_type(-1);
  bool hardEmission=false;
  if(wgtb>UseRandom::rnd()*(weights_[0]+wgtb)) {
    wgtb *= UseRandom::rnd();
    unsigned int itype=1;
    for(;itype<weights_.size();++itype) {
      if(weights_[itype]>=wgtb) break;
      wgtb-=weights_[itype];
    }
    emission_type = itype;
    hardEmission = true;
  }
  // generate a hard emission from the sudakov
  else {
    int qcd_type(-1),qed_type(-1);
    Energy qcd_pT(-GeV),qed_pT(-GeV);
    for(unsigned int ix=0;ix<inter.size();++ix) {
      if(inter[ix]==ShowerInteraction::QCD)
	hardQCDEmission(particlesToShower,qcd_type,qcd_pT);
      else if(inter[ix]==ShowerInteraction::QED)
	hardQEDEmission(particlesToShower,qed_type,qed_pT);
    }
    if(qcd_type<0&&qed_type<0) {
      Energy ptminQED= min(min(minpTQEDII_,minpTQEDIF_),minpTQEDFF_);
      for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
	particlesToShower[ix]->maximumpT(minpTQCD_,ShowerInteraction::QCD);
	particlesToShower[ix]->maximumpT(ptminQED ,ShowerInteraction::QED);
      }
      return HardTreePtr();
    }
    else if(qcd_type>0&&qed_type>0) {
      if(qcd_pT>=qed_pT) {
	emission_type = qcd_type;
      }
      else {
	emission_type = qed_type;
      }
    }
    else if(qcd_type>0) {
      emission_type = qcd_type;
    }
    else {
      emission_type = qed_type;
    }
  }
  // construct the HardTree object needed to perform the showers
  ShowerParticleVector newparticles;
  // make the particles for the HardTree
  // create the partons
  // q qbar -> g X
  vector<Lorentz5Momentum> pnew;
  if(emission_type==1||emission_type==3) {
    newparticles.push_back(new_ptr(ShowerParticle(_partons[0]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[1]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(gluon_           , true)));
    if(emission_type ==1) pnew = realEmissionQCDGluon1_;
    else                  pnew = realEmissionQCDGluon2_;
  }
  // q g    -> q X
  else if(emission_type==4) {
    newparticles.push_back(new_ptr(ShowerParticle(_partons[0]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(gluon_           ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[1]->CC(), true)));
    pnew = realEmissionQCDQuark2_;
  }
  // g qbar -> qbar X
  else if(emission_type==2) {
    newparticles.push_back(new_ptr(ShowerParticle(gluon_           ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[1]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[0]->CC(), true)));
    pnew = realEmissionQCDQuark1_;
  }
  else if(emission_type==5||emission_type==7) {
    newparticles.push_back(new_ptr(ShowerParticle(_partons[0]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[1]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(gamma_           , true)));
    if(emission_type ==5) pnew = realEmissionQEDPhoton1_;
    else                  pnew = realEmissionQEDPhoton2_;
  }
  else if(emission_type==8) {
    newparticles.push_back(new_ptr(ShowerParticle(_partons[0]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(gamma_           ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[1]->CC(), true)));
    pnew = realEmissionQEDQuark2_;
  }
  else if(emission_type==6) {
    newparticles.push_back(new_ptr(ShowerParticle(gamma_           ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[1]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[0]->CC(), true)));
    pnew = realEmissionQEDQuark1_;
  }
  else if(emission_type==9 || emission_type==10) {
    newparticles.push_back(new_ptr(ShowerParticle(_partons[0]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[1]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(gamma_           , true)));
    if(emission_type ==9) pnew = realEmissionQEDFFLepton3_;
    else                  pnew = realEmissionQEDFFLepton4_;
  }
  else if(emission_type>10&&emission_type<=18) {
    newparticles.push_back(new_ptr(ShowerParticle(_partons[0]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[1]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(gamma_           , true)));
    if     (emission_type==11||emission_type==15) pnew =realEmissionQEDIF13_;
    else if(emission_type==12||emission_type==16) pnew =realEmissionQEDIF14_;
    else if(emission_type==13||emission_type==17) pnew =realEmissionQEDIF23_;
    else if(emission_type==14||emission_type==18) pnew =realEmissionQEDIF24_;
  }
  else {
    assert(false);
  }
  // transform to get directions right if hard emission
  LorentzRotation R;
  if(!_quarkplus&&hardEmission) {
    PPair partons = make_pair(particlesToShower[0]->progenitor(),
			      particlesToShower[1]->progenitor());
    if(partons.first->id()<0) swap(partons.first,partons.second);
    Boost bv = (partons.first->momentum()+
		partons.second->momentum()).boostVector();
    R = LorentzRotation(-bv);
    R.rotateY(-partons.first->momentum().theta());
    R.boost(bv);
    for(unsigned int ix=0;ix<pnew.size();++ix)
      pnew[ix].transform(R);
  }
  // work out which particle emits
  int iemit(-1),ispect(-1);
  if(emission_type==1 || emission_type==2 || 
     emission_type==5 || emission_type==6) {
    iemit  = 0;
    ispect = 1;
  }
  else if(emission_type==3 || emission_type==4 || 
	  emission_type==7 || emission_type==8) {
    iemit  = 1;
    ispect = 0;
  }
  else if(emission_type==9) {
    iemit  = 3;
    ispect = 4;
  }
  else if(emission_type==10) {
    iemit  = 4;
    ispect = 3;
  }
  else if(emission_type>10&&emission_type<=18) {
    iemit  = (emission_type==11||emission_type==15||
	      emission_type==12||emission_type==16) ? 0 : 1;
    ispect = (emission_type==11||emission_type==15||
	      emission_type==13||emission_type==17) ? 3 : 4;
    if(emission_type>14) swap(iemit,ispect);
  }
  assert(iemit>=0&&ispect>=0);
  // work out which interaction
  ShowerInteraction::Type interaction = 
    emission_type<=4 ? ShowerInteraction::QCD : ShowerInteraction::QED;
  // set the momenta
  for(unsigned int ix=0;ix<2;++ix) newparticles[ix]->set5Momentum(pnew[ix]);
  newparticles[2]->set5Momentum(pnew[4]);
  // ISR
  if(iemit<2) {
    // create the off-shell particle
    Lorentz5Momentum poff=pnew[iemit]-pnew[4];
    poff.rescaleMass();
    newparticles.push_back(new_ptr(ShowerParticle(_partons[iemit],false)));
    newparticles.back()->set5Momentum(poff);
  }
  // FSR
  else {
    // create the off-shell particle
    Lorentz5Momentum poff=pnew[iemit-1]+pnew[4];
    poff.rescaleMass();
    newparticles.push_back(new_ptr(ShowerParticle(mePartonData()[iemit-1],true)));
    newparticles.back()->set5Momentum(poff);
  }
  // outgoing particles
  for(unsigned int ix=2;ix<4;++ix) {
    newparticles.push_back(new_ptr(ShowerParticle(particlesToShower[ix]
						  ->progenitor()->dataPtr(),
						  true)));
  }
  newparticles[4]->set5Momentum(pnew[2]);
  newparticles[5]->set5Momentum(pnew[3]);
  vector<HardBranchingPtr> inBranch,hardBranch;
  // create the branchings for the incoming particles
  inBranch.push_back(new_ptr(HardBranching(newparticles[0],SudakovPtr(),
					   HardBranchingPtr(),HardBranching::Incoming)));
  inBranch.push_back(new_ptr(HardBranching(newparticles[1],SudakovPtr(),
					   HardBranchingPtr(),HardBranching::Incoming)));
  // ISR
  if(iemit<2) {
    // intermediate IS particle
    hardBranch.push_back(new_ptr(HardBranching(newparticles[3],SudakovPtr(),
					       inBranch[iemit],HardBranching::Incoming)));
    inBranch[iemit]->addChild(hardBranch.back());
    // create the branching for the emitted jet
    inBranch[iemit]->addChild(new_ptr(HardBranching(newparticles[2],SudakovPtr(),
						    inBranch[iemit],
						    HardBranching::Outgoing)));
    // add other particle
    hardBranch.push_back(inBranch[iemit==0 ? 1 : 0]);
    // outgoing particles
    for(unsigned int ix=4;ix<newparticles.size();++ix) {
      hardBranch.push_back(new_ptr(HardBranching(newparticles[ix],SudakovPtr(),
						 HardBranchingPtr(),
						 HardBranching::Outgoing)));
    }
  }
  // FSR
  else {
    hardBranch.push_back(inBranch[0]);
    hardBranch.push_back(inBranch[1]);
    if(iemit==3) {
      hardBranch.push_back(new_ptr(HardBranching(newparticles[3],SudakovPtr(),
						 HardBranchingPtr(),
						 HardBranching::Outgoing)));
      hardBranch.back()->addChild(new_ptr(HardBranching(newparticles[4],SudakovPtr(),
							hardBranch.back(),
							HardBranching::Outgoing)));
      hardBranch.back()->addChild(new_ptr(HardBranching(newparticles[2],SudakovPtr(),
							hardBranch.back(),
							HardBranching::Outgoing)));
      hardBranch.push_back(new_ptr(HardBranching(newparticles[5],SudakovPtr(),
						 HardBranchingPtr(),
						 HardBranching::Outgoing)));
	
    }
    else {
      hardBranch.push_back(new_ptr(HardBranching(newparticles[4],SudakovPtr(),
						 HardBranchingPtr(),
						 HardBranching::Outgoing)));
      hardBranch.push_back(new_ptr(HardBranching(newparticles[3],SudakovPtr(),
						 HardBranchingPtr(),
						 HardBranching::Outgoing)));
      hardBranch.back()->addChild(new_ptr(HardBranching(newparticles[5],SudakovPtr(),
							hardBranch.back(),
							HardBranching::Outgoing)));
      hardBranch.back()->addChild(new_ptr(HardBranching(newparticles[2],SudakovPtr(),
							hardBranch.back(),
							HardBranching::Outgoing)));
    }
  }
  // set the evolution partners
  if(emission_type<=4) {
    hardBranch[0]->colourPartner(hardBranch[1]);
    hardBranch[1]->colourPartner(hardBranch[0]);
  }
  else if(emission_type<=10) {
    hardBranch[2]->colourPartner(hardBranch[3]);
    hardBranch[3]->colourPartner(hardBranch[2]);
    hardBranch[2]->colourPartner(hardBranch[3]);
    hardBranch[3]->colourPartner(hardBranch[2]);
  }
  else if(emission_type<=18) {
    if(iemit>2) {
      hardBranch[iemit -1]->colourPartner(hardBranch[ispect  ]);
      hardBranch[ispect  ]->colourPartner(hardBranch[iemit -1]);
    }
    else {
      hardBranch[iemit   ]->colourPartner(hardBranch[ispect-1]);
      hardBranch[ispect-1]->colourPartner(hardBranch[iemit   ]);
    }
  }
  else {
    assert(false);
  }
  // make the tree
  HardTreePtr hardtree=new_ptr(HardTree(hardBranch,inBranch,interaction));
  // connect the ShowerParticles with the branchings
  // and set the maximum pt for the radiation
  Energy pt = pnew[4].perp();
  set<HardBranchingPtr> hard=hardtree->branchings();
  Energy ptminQED= min(min(minpTQEDII_,minpTQEDIF_),minpTQEDFF_);
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    particlesToShower[ix]->maximumpT(max(pt,minpTQCD_),ShowerInteraction::QCD);
    particlesToShower[ix]->maximumpT(max(pt,ptminQED),ShowerInteraction::QED);
    for(set<HardBranchingPtr>::const_iterator mit=hard.begin();
	mit!=hard.end();++mit) {
      if(particlesToShower[ix]->progenitor()->id()==(*mit)->branchingParticle()->id()&&
	 (( particlesToShower[ix]->progenitor()->isFinalState()&&
	    (**mit).status()==HardBranching::Outgoing)||
	  (!particlesToShower[ix]->progenitor()->isFinalState()&&
	   (**mit).status()==HardBranching::Incoming))) {
	hardtree->connect(particlesToShower[ix]->progenitor(),*mit);
	if((**mit).status()==HardBranching::Incoming) {
	  (*mit)->beam(particlesToShower[ix]->original()->parents()[0]);
	  HardBranchingPtr parent=(*mit)->parent();
	  while(parent) {
	    parent->beam(particlesToShower[ix]->original()->parents()[0]);
	    parent=parent->parent();
	  };
	}
      }
    }
  }
  ColinePtr newline=new_ptr(ColourLine());
  for(set<HardBranchingPtr>::const_iterator cit=hardtree->branchings().begin();
      cit!=hardtree->branchings().end();++cit) {
    if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3)
      newline->addColoured((**cit).branchingParticle());
    else if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3bar)
      newline->addAntiColoured((**cit).branchingParticle());
  }
  // return the tree
//   generator()->log() << "Had hard radiation " << iemit << "\n";
//   generator()->log() << *hardtree << "\n";
//   cerr << "Had hard radiation " << iemit << "\n";
//   cerr << *hardtree << "\n";
  return hardtree;
}
