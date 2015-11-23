// -*- C++ -*-
//
// MPIHandler.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MPIHandler class.
//

#include "MPIHandler.h"

#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Handlers/SubProcessHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/Handlers/CascadeHandler.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Cuts/JetCuts.h"
#include "ThePEG/Cuts/SimpleKTCut.h"
#include "ThePEG/PDF/PartonExtractor.h"

#include "gsl/gsl_sf_bessel.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MPIHandler * MPIHandler::currentHandler_ = 0;

bool MPIHandler::beamOK() const {
  return (HadronMatcher::Check(*eventHandler()->incoming().first)  &&
	  HadronMatcher::Check(*eventHandler()->incoming().second) );
}


tStdXCombPtr MPIHandler::generate(unsigned int sel) {
  //generate a certain process
  if(sel+1 > processHandlers().size())
    throw Exception() << "MPIHandler::generate called with argument out of range"
                      << Exception::runerror;

  return processHandlers()[sel]->generate();
}

IBPtr MPIHandler::clone() const {
  return new_ptr(*this);
}

IBPtr MPIHandler::fullclone() const {
  return new_ptr(*this);
}

void MPIHandler::finalize() {
  if( beamOK() ){
    statistics();
  }
}

void MPIHandler::initialize() {
  currentHandler_ = this;
  useMe();
  theHandler = generator()->currentEventHandler(); 
  //stop if the EventHandler is not present:
  assert(theHandler);

  //check if MPI is wanted
  if( !beamOK() ){
    throw Exception()  << "You have requested multiple parton-parton scattering,\n"
		       << "but the model is not forseen for the beam setup you chose.\n" 
		       << "You should therefore disable that by setting XXXGenerator:EventHandler:"
		       << "CascadeHandler:MPIHandler to NULL"
                       << Exception::runerror;
  }

  numSubProcs_ = subProcesses().size();

  if( numSubProcs_ != cuts().size() ) 
    throw Exception() << "MPIHandler::each SubProcess needs a Cuts Object"
		      << "ReferenceVectors are not equal in size"
		      << Exception::runerror;

  if( additionalMultiplicities_.size()+1 != numSubProcs_ )
    throw Exception() << "MPIHandler: each additional SubProcess needs "
		      << "a multiplicity assigned. This can be done in with "
		      << "insert MPIHandler:additionalMultiplicities 0 1"
		      << Exception::runerror;

  //identicalToUE_ = 0 hard process is identical to ue, -1 no one
  if( identicalToUE_ > (int)numSubProcs_ || identicalToUE_ < -1 )
    throw Exception() << "MPIHandler:identicalToUE has disallowed value"
		      << Exception::runerror;

  // override the cuts for the additional scatters if energyExtrapolation_ is
  // set
  if (energyExtrapolation_ != 0 ) {
    overrideUECuts();
  }

  tcPDPtr gluon=getParticleData(ParticleID::g);
  //determine ptmin
  Ptmin_ = cuts()[0]->minKT(gluon);

  if(identicalToUE_ == -1){
    algorithm_ = 2;
  }else{
    if(identicalToUE_ == 0){
      //Need to work a bit, in case of LesHouches events for QCD2to2
      Ptr<StandardEventHandler>::pointer eH =
	dynamic_ptr_cast<Ptr<StandardEventHandler>::pointer>(eventHandler());

      if( eH ) {
	PtOfQCDProc_ = eH->cuts()->minKT(gluon);
	// find the jet cut in the new style cuts
	for(unsigned int ix=0;ix<eH->cuts()->multiCuts().size();++ix) {
	  Ptr<JetCuts>::pointer jetCuts =
	    dynamic_ptr_cast<Ptr<JetCuts>::pointer>(eH->cuts()->multiCuts()[ix]);
	  if(jetCuts) {
	    Energy ptMin=1e30*GeV;
	    for(unsigned int iy=0;iy<jetCuts->jetRegions().size();++iy) {
	      ptMin = min(ptMin,jetCuts->jetRegions()[iy]->ptMin());
	    }
	    if(ptMin<1e29*GeV&&ptMin>PtOfQCDProc_) PtOfQCDProc_ = ptMin;
	  }
	}
      }
      else {
	if(PtOfQCDProc_ == -1.0*GeV)
	  throw Exception() << "MPIHandler: You need to specify the pt cutoff "
			    << "used to in the LesHouches file for QCD2to2 events"
			    << Exception::runerror;
      }
    }
    else {
      PtOfQCDProc_ = cuts()[identicalToUE_]->minKT(gluon);
    }

    if(PtOfQCDProc_ > 2*Ptmin_)
      algorithm_ = 1;
    else
      algorithm_ = 0;

    if(PtOfQCDProc_ == ZERO)//pure MinBias mode
      algorithm_ = -1;
  }

  //Init all subprocesses
  for(unsigned int i=0; i<numSubProcs_; i++){
    theProcessHandlers.push_back(new_ptr(ProcessHandler()));
    processHandlers().back()->initialize(subProcesses()[i], 
					 cuts()[i], eventHandler());
    processHandlers().back()->initrun();
  }

  //now calculate the individual Probabilities
  XSVector UEXSecs;
  UEXSecs.push_back(processHandlers()[0]->integratedXSec());
  //save the hard cross section
  hardXSec_ = UEXSecs.front();

  //determine sigma_soft and beta
  if(softInt_){//check that soft ints are requested
    GSLBisection rootFinder;

    if(twoComp_){
      //two component model
      /*
      GSLMultiRoot eqSolver;
      slopeAndTotalXSec eq(this);
      pair<CrossSection, Energy2> res = eqSolver.value(eq, 10*millibarn, 0.6*GeV2);
      softXSec_ = res.first;
      softMu2_ = res.second;
      */
      slopeBisection fs(this);
      try{
	softMu2_ = rootFinder.value(fs, 0.3*GeV2, 1.*GeV2);
	softXSec_ = fs.softXSec();
      }catch(GSLBisection::IntervalError){
	throw Exception() <<
        "\n**********************************************************\n"
	  "* Inconsistent MPI parameter choice for this beam energy *\n"
	  "**********************************************************\n"
	  "MPIHandler parameter choice is unable to reproduce\n"
	  "the total cross section. Please check arXiv:0806.2949\n"
	  "for the allowed parameter space."
			  << Exception::runerror;
      }

    }else{
      //single component model
      TotalXSecBisection fn(this);
      try{
	softXSec_ = rootFinder.value(fn, 0*millibarn, 5000*millibarn);
      }catch(GSLBisection::IntervalError){
	throw Exception() <<
	"\n**********************************************************\n"
	  "* Inconsistent MPI parameter choice for this beam energy *\n"
	  "**********************************************************\n"
	  "MPIHandler parameter choice is unable to reproduce\n"
	  "the total cross section. Please check arXiv:0806.2949\n"
	  "for the allowed parameter space."
			  << Exception::runerror;      
      }
    }

    //now get the differential cross section at ptmin
    ProHdlPtr qcd = new_ptr(ProcessHandler());
    Energy eps = 0.1*GeV;
    Energy ptminPlus = Ptmin_ + eps;
    
    Ptr<SimpleKTCut>::pointer ktCut = new_ptr(SimpleKTCut(ptminPlus));
    ktCut->init();
    ktCut->initrun();

    CutsPtr qcdCut = new_ptr(Cuts(2*ptminPlus));
    qcdCut->add(dynamic_ptr_cast<tOneCutPtr>(ktCut));
    qcdCut->init();
    qcdCut->initrun();
      
    qcd->initialize(subProcesses()[0], qcdCut, eventHandler());
    qcd->initrun();
    
    // ds/dp_T^2 = 1/2/p_T ds/dp_T
    DiffXSec hardPlus = (hardXSec_-qcd->integratedXSec())/(2*Ptmin_*eps);
    
    betaBisection fn2(softXSec_, hardPlus, Ptmin_);
    try{
      beta_ = rootFinder.value(fn2, -10/GeV2, 2/GeV2);
    }catch(GSLBisection::IntervalError){
      throw Exception() << "MPIHandler: slope of soft pt spectrum couldn't be "
			<< "determined."
			<< Exception::runerror;    
    }
  }

  Probs(UEXSecs);
  //MultDistribution("probs.test");

  UEXSecs.clear();
}

void MPIHandler::MultDistribution(string filename) const {
  ofstream file;
  double p(0.0), pold(0.0);
  file.open(filename.c_str());
  //theMultiplicities  
  Selector<MPair>::const_iterator it = theMultiplicities.begin();

  while(it != theMultiplicities.end()){
    p = it->first;
    file << it->second.first << " " << it->second.second
	 << " " << p-pold << '\n';
    it++;
    pold = p;
  }

  file << "sum of all probabilities: " << theMultiplicities.sum() 
       << endl;

  file.close();
}

void MPIHandler::statistics() const {
  ostream & file = generator()->misc();
  
  string line = "======================================="
    "=======================================\n";

  for(unsigned int i=0; i<cuts().size(); i++){
    Stat tot;
    if(i == 0)
      file << "Statistics for the UE process: \n";
    else
      file << "Statistics for additional hard Process " << i << ": \n";

    processHandlers()[i]->statistics(file, tot);
    file << "\n";
  }

  if(softInt_){
    file << line
	 << "Eikonalized and soft cross sections:\n\n"
	 << "Model parameters:                    "
	 << "ptmin:   " << Ptmin_/GeV << " GeV"
      	 << ", mu2: " << invRadius_/sqr(1.*GeV) << " GeV2\n"
	 << "                                     "
	 << "DL mode: " << DLmode_
	 << ", CMenergy: " << generator()->maximumCMEnergy()/GeV
	 << " GeV" << '\n'
	 << "hard inclusive cross section (mb):   "
	 << hardXSec_/millibarn << '\n'
	 << "soft inclusive cross section (mb):   "
	 << softXSec_/millibarn << '\n'
	 << "total cross section (mb):            "
	 << totalXSecExp()/millibarn << '\n'
	 << "inelastic cross section (mb):        "
	 << inelXSec_/millibarn << '\n'
	 << "soft inv radius (GeV2):              "
	 << softMu2_/GeV2 << '\n'
	 << "slope of soft pt spectrum (1/GeV2):  "
	 << beta_*sqr(1.*GeV) << '\n'
	 << "Average hard multiplicity:           "
	 << avgNhard_ << '\n'
	 << "Average soft multiplicity:           "
	 << avgNsoft_ << '\n' << line << endl;
  }else{
    file << line
	 << "Eikonalized and soft cross sections:\n\n"
	 << "Model parameters:                    "
	 << "ptmin:   " << Ptmin_/GeV << " GeV"
      	 << ", mu2: " << invRadius_/sqr(1.*GeV) << " GeV2\n"
	 << "                                     "
	 << ", CMenergy: " << generator()->maximumCMEnergy()/GeV
	 << " GeV" << '\n'
	 << "hard inclusive cross section (mb):   "
	 << hardXSec_/millibarn << '\n'
	 << "Average hard multiplicity:           "
	 << avgNhard_ << '\n' << line << endl;
  }
}

unsigned int MPIHandler::multiplicity(unsigned int sel){
  if(sel==0){//draw from the pretabulated distribution
    MPair m = theMultiplicities.select(UseRandom::rnd());
    softMult_ = m.second;
    return m.first;
  } else{ //fixed multiplicities for the additional hard scatters
    if(additionalMultiplicities_.size() < sel)
      throw Exception() << "MPIHandler::multiplicity: process index "
			<< "is out of range"
			<< Exception::runerror;

    return additionalMultiplicities_[sel-1];
  }
}


void MPIHandler::Probs(XSVector UEXSecs) {
  GSLIntegrator integrator;
  unsigned int iH(1), iS(0);
  double P(0.0);
  double P0(0.0);//the probability for i hard and zero soft scatters
  Length bmax(500.0*sqrt(millibarn));
  //only one UE process will be drawn from a probability distribution,
  //so check that.
  assert(UEXSecs.size() == 1);

  for ( XSVector::const_iterator it = UEXSecs.begin();
        it != UEXSecs.end(); ++it ) {

    iH = 0; 

    //get the inel xsec
    Eikonalization inelint(this, -1, *it, softXSec_, softMu2_);
    inelXSec_ = integrator.value(inelint, ZERO, bmax);

    avgNhard_ = 0.0;
    avgNsoft_ = 0.0;
    bmax = 10.0*sqrt(millibarn);
    do{//loop over hard ints
      if(Algorithm()>-1 && iH==0){
	iH++;
	continue;
      }
      iS = 0;
      do{//loop over soft ints
	if( ( Algorithm() == -1 && iS==0 && iH==0 ) ){
	  iS++;
	  continue;
	}

	Eikonalization integrand(this, iH*100+iS, *it, softXSec_, softMu2_);
      
	if(Algorithm() > 0){
	  P = integrator.value(integrand, ZERO, bmax) / *it;
	}else{
	  P = integrator.value(integrand, ZERO, bmax) / inelXSec_;
	}
	//store the probability
	if(Algorithm()>-1){
	  theMultiplicities.insert(P, make_pair(iH-1, iS));
	  avgNhard_ += P*(iH-1);
	}else{
	  theMultiplicities.insert(P, make_pair(iH, iS));
	  avgNhard_ += P*(iH);
	}
	avgNsoft_ += P*iS;
	if(iS==0)
	  P0 = P;

	iS++;
      } while ( (iS < maxScatters_) && (iS < 5 || P > 1.e-15 ) && softInt_ );
      iH++;
    } while ( (iH < maxScatters_) && (iH < 5 || P0 > 1.e-15) );
  }
}


// calculate the integrand
Length Eikonalization::operator() (Length b) const {
  unsigned int Nhard(0), Nsoft(0);
  //fac is just: d^2b=fac*db despite that large number
  Length fac(Constants::twopi*b);
  
  double chiTot(( theHandler->OverlapFunction(b)*hardXSec_ + 
		  theHandler->OverlapFunction(b, softMu2_)*softXSec_ ) / 2.0);

  //total cross section wanted
  if(theoption == -2) return 2 * fac * ( 1 - exp(-chiTot) );
  //inelastic cross section
  if(theoption == -1) return fac * ( 1 - exp(- 2.0 * chiTot) );

  if(theoption >= 0){
    //encode multiplicities as: N_hard*100 + N_soft   
    Nhard = theoption/100;
    Nsoft = theoption%100;

    if(theHandler->Algorithm() > 0){
      //P_n*sigma_hard: n-1 extra scatters + 1 hard scatterer != hardXSec_
      return fac * Nhard * theHandler->poisson(b, hardXSec_, Nhard) *
	theHandler->poisson(b, softXSec_, Nsoft, softMu2_);
    }else{
      //P_n*sigma_inel: n extra scatters
      return fac * theHandler->poisson(b, hardXSec_, Nhard) *
	theHandler->poisson(b, softXSec_, Nsoft, softMu2_);
    }

  }else{
    throw Exception() << "Parameter theoption in Struct Eikonalization in " 
		      << "MPIHandler.cc has not allowed value"
                      << Exception::runerror;
    return 0.0*meter;
  }
}

InvEnergy2 slopeBisection::operator() (Energy2 softMu2) const {
  GSLBisection root;
  //single component model
  TotalXSecBisection fn(handler_, softMu2);
  try{
    softXSec_ = root.value(fn, 0*millibarn, 5000*millibarn);
  }catch(GSLBisection::IntervalError){
    throw Exception() << "MPIHandler 2-Component model didn't work out."
		      << Exception::runerror;
  }
  return handler_->slopeDiff(softXSec_, softMu2);
}

LengthDiff slopeInt::operator() (Length b) const {
  //fac is just: d^2b=fac*db
  Length fac(Constants::twopi*b);
  
  double chiTot(( handler_->OverlapFunction(b)*hardXSec_ + 
		  handler_->OverlapFunction(b, softMu2_)*softXSec_ ) / 2.0);

  InvEnergy2 b2 = sqr(b/hbarc);
  //B*sigma_tot
  return fac * b2 * ( 1 - exp(-chiTot) );
}

double MPIHandler::factorial (unsigned int n) const {

  double f[] = {1.,1.,2.,6.,24.,120.,720.,5040.,40320.,362880.,3.6288e6,
		3.99168e7,4.790016e8,6.2270208e9,8.71782912e10,1.307674368e12,
		2.0922789888e13,3.55687428096e14,6.402373705728e15,1.21645100408832e17,
		2.43290200817664e18,5.10909421717094e19,1.12400072777761e21,
		2.5852016738885e22,6.20448401733239e23,1.5511210043331e25,
		4.03291461126606e26,1.08888694504184e28,3.04888344611714e29,
		8.8417619937397e30,2.65252859812191e32,8.22283865417792e33,
		2.63130836933694e35,8.68331761881189e36,2.95232799039604e38,
		1.03331479663861e40,3.71993326789901e41,1.37637530912263e43,
		5.23022617466601e44,2.03978820811974e46,8.15915283247898e47,
		3.34525266131638e49,1.40500611775288e51,6.04152630633738e52,
		2.65827157478845e54,1.1962222086548e56,5.50262215981209e57,
		2.58623241511168e59,1.24139155925361e61,6.08281864034268e62,
		3.04140932017134e64,1.55111875328738e66,8.06581751709439e67,
		4.27488328406003e69,2.30843697339241e71,1.26964033536583e73,
		7.10998587804863e74,4.05269195048772e76,2.35056133128288e78,
		1.3868311854569e80,8.32098711274139e81,5.07580213877225e83,
		3.14699732603879e85,1.98260831540444e87,1.26886932185884e89,
		8.24765059208247e90,5.44344939077443e92,3.64711109181887e94,
		2.48003554243683e96,1.71122452428141e98,1.19785716699699e100,
		8.50478588567862e101,6.12344583768861e103,4.47011546151268e105,
		3.30788544151939e107,2.48091408113954e109,1.88549470166605e111,
		1.45183092028286e113,1.13242811782063e115,8.94618213078298e116,
		7.15694570462638e118,5.79712602074737e120,4.75364333701284e122,
		3.94552396972066e124,3.31424013456535e126,2.81710411438055e128,
		2.42270953836727e130,2.10775729837953e132,1.85482642257398e134,
		1.65079551609085e136,1.48571596448176e138,1.3520015276784e140,
		1.24384140546413e142,1.15677250708164e144,1.08736615665674e146,
		1.03299784882391e148,9.9167793487095e149,9.61927596824821e151,
		9.42689044888325e153,9.33262154439442e155,9.33262154439442e157};

  if(n > maxScatters_) 
        throw Exception() << "MPIHandler::factorial called with too large argument"
                      << Exception::runerror;
  else
    return f[n];
}

InvArea MPIHandler::OverlapFunction(Length b, Energy2 mu2) const {
  if(mu2 == ZERO)
    mu2 = invRadius_;

  InvLength mu = sqrt(mu2)/hbarc;
  return (sqr(mu)/96/Constants::pi)*pow(mu*b, 3)*(gsl_sf_bessel_Kn(3, mu*b));
}

double MPIHandler::poisson(Length b, CrossSection sigma, unsigned int N, Energy2 mu2) const {
  if(sigma > 0*millibarn){
    return pow(OverlapFunction(b, mu2)*sigma, (double)N)/factorial(N)
      *exp(-OverlapFunction(b, mu2)*sigma);
  }else{
    return (N==0) ? 1.0 : 0.0;
  }
}

CrossSection MPIHandler::totalXSecDiff(CrossSection softXSec, 
				       Energy2 softMu2) const {
  GSLIntegrator integrator;
  Eikonalization integrand(this, -2, hardXSec_, softXSec, softMu2);
  Length bmax = 500.0*sqrt(millibarn);

  CrossSection tot = integrator.value(integrand, ZERO, bmax);
  return (tot-totalXSecExp());
}

InvEnergy2 MPIHandler::slopeDiff(CrossSection softXSec, 
				 Energy2 softMu2) const {
  GSLIntegrator integrator;
  Eikonalization integrand(this, -2, hardXSec_, softXSec, softMu2);
  Length bmax = 500.0*sqrt(millibarn);

  CrossSection tot = integrator.value(integrand, ZERO, bmax);
  
  slopeInt integrand2(this, hardXSec_, softXSec, softMu2);
  
  return integrator.value(integrand2, ZERO, bmax)/tot - slopeExp();
}

CrossSection MPIHandler::totalXSecExp() const {
  if(totalXSecExp_ != 0*millibarn)
    return totalXSecExp_;

  double pom_old = 0.0808;
  CrossSection coef_old = 21.7*millibarn;
  double pom_new_hard = 0.452;
  CrossSection coef_new_hard = 0.0139*millibarn;
  double pom_new_soft = 0.0667;
  CrossSection coef_new_soft = 24.22*millibarn;

  Energy energy(generator()->maximumCMEnergy());

  switch(DLmode_){
  case 1://old DL extrapolation
    return coef_old * pow(energy/GeV, 2*pom_old);
    break;

  case 2://old DL extrapolation fixed to CDF
    return 81.8*millibarn * pow(energy/1800.0/GeV, 2*pom_old);
    break;
    
  case 3://new DL extrapolation
    return coef_new_hard * pow(energy/GeV, 2*pom_new_hard) + 
      coef_new_soft * pow(energy/GeV, 2*pom_new_soft);
    break;
    
  default:
    throw Exception() << "MPIHandler::totalXSecExp non-existing mode selected"
                      << Exception::runerror;   
  }
}

InvEnergy2 MPIHandler::slopeExp() const{
  //Currently return the slope as calculated in the pomeron fit by
  //Donnachie & Landshoff
  Energy energy(generator()->maximumCMEnergy());
  //slope at
  Energy e_0 = 1800*GeV;
  InvEnergy2 b_0 = 17/GeV2;
  return b_0 + log(energy/e_0)/GeV2;
}

void MPIHandler::overrideUECuts() {
  if(energyExtrapolation_==1)
    Ptmin_ = EEparamA_ * log(generator()->maximumCMEnergy() / EEparamB_);
  else if(energyExtrapolation_==2)
    Ptmin_ = pT0_*pow(double(generator()->maximumCMEnergy()/refScale_),b_);
  else
    assert(false);
  // create a new SimpleKTCut object with the calculated ptmin value
  Ptr<SimpleKTCut>::pointer newUEktCut = new_ptr(SimpleKTCut(Ptmin_));
  newUEktCut->init();
  newUEktCut->initrun();

  // create a new Cuts object with MHatMin = 2 * Ptmin_
  CutsPtr newUEcuts = new_ptr(Cuts(2*Ptmin_));
  newUEcuts->add(dynamic_ptr_cast<tOneCutPtr>(newUEktCut));
  newUEcuts->init();
  newUEcuts->initrun();

  // replace the old Cuts object
  cuts()[0] = newUEcuts;
}

void MPIHandler::persistentOutput(PersistentOStream & os) const {
  os << theMultiplicities << theHandler
     << theSubProcesses << theCuts << theProcessHandlers
     << additionalMultiplicities_ << identicalToUE_ 
     << ounit(PtOfQCDProc_, GeV) << ounit(Ptmin_, GeV) 
     << ounit(hardXSec_, millibarn) << ounit(softXSec_, millibarn)
     << ounit(beta_, 1/GeV2)
     << algorithm_ << ounit(invRadius_, GeV2)
     << numSubProcs_ << colourDisrupt_ << softInt_ << twoComp_ 
     << DLmode_ << ounit(totalXSecExp_, millibarn)
     << energyExtrapolation_ << ounit(EEparamA_, GeV) << ounit(EEparamB_, GeV)
     << ounit(refScale_,GeV) << ounit(pT0_,GeV) << b_
     << avgNhard_ << avgNsoft_ << softMult_ 
     << ounit(inelXSec_, millibarn) 
     << ounit(softMu2_, GeV2);
}

void MPIHandler::persistentInput(PersistentIStream & is, int) {
  is >> theMultiplicities >> theHandler
     >> theSubProcesses >> theCuts >> theProcessHandlers
     >> additionalMultiplicities_ >> identicalToUE_ 
     >> iunit(PtOfQCDProc_, GeV) >> iunit(Ptmin_, GeV)
     >> iunit(hardXSec_, millibarn) >> iunit(softXSec_, millibarn)
     >> iunit(beta_, 1/GeV2)
     >> algorithm_ >> iunit(invRadius_, GeV2)
     >> numSubProcs_ >> colourDisrupt_ >> softInt_ >> twoComp_ 
     >> DLmode_ >> iunit(totalXSecExp_, millibarn)
     >> energyExtrapolation_ >> iunit(EEparamA_, GeV) >> iunit(EEparamB_, GeV)
     >> iunit(refScale_,GeV) >> iunit(pT0_,GeV) >> b_
     >> avgNhard_ >> avgNsoft_ >> softMult_ 
     >> iunit(inelXSec_, millibarn) 
     >> iunit(softMu2_, GeV2);
  currentHandler_ = this;
}

void MPIHandler::clean() {
  // ThePEG's event handler's usual event cleanup doesn't reach these
  // XCombs. Need to do it by hand here.
  for ( size_t i = 0; i < theSubProcesses.size(); ++i ) {
    theSubProcesses[i]->pExtractor()->lastXCombPtr()->clean();
  }
}


ClassDescription<MPIHandler> MPIHandler::initMPIHandler;
// Definition of the static class description member.

void MPIHandler::Init() {

  static ClassDocumentation<MPIHandler> documentation
    ("The MPIHandler class is the main administrator of the multiple interaction model", 
     "The underlying event was simulated with an eikonal model for multiple partonic interactions."
     "Details can be found in Ref.~\\cite{Bahr:2008dy,Bahr:2009ek}.", 
     "%\\cite{Bahr:2008dy}\n"
     "\\bibitem{Bahr:2008dy}\n"
     "  M.~Bahr, S.~Gieseke and M.~H.~Seymour,\n"
     "  ``Simulation of multiple partonic interactions in Herwig,''\n"
     "  JHEP {\\bf 0807}, 076 (2008)\n"
     "  [arXiv:0803.3633 [hep-ph]].\n"
     "  %%CITATION = JHEPA,0807,076;%%\n"
    "\\bibitem{Bahr:2009ek}\n"
     "  M.~Bahr, J.~M.~Butterworth, S.~Gieseke and M.~H.~Seymour,\n"
     "  ``Soft interactions in Herwig,''\n"
     "  arXiv:0905.4671 [hep-ph].\n"
     "  %%CITATION = ARXIV:0905.4671;%%\n"
     );
  
  static RefVector<MPIHandler,SubProcessHandler> interfaceSubhandlers
    ("SubProcessHandlers",
     "The list of sub-process handlers used in this EventHandler. ",
     &MPIHandler::theSubProcesses, -1, false, false, true, false, false);

  static RefVector<MPIHandler,Cuts> interfaceCuts
    ("Cuts",
     "List of cuts used for the corresponding list of subprocesses. These cuts "
     "should not be overidden in individual sub-process handlers.",
     &MPIHandler::theCuts, -1, false, false, true, false, false);

  static Parameter<MPIHandler,Energy2> interfaceInvRadius
    ("InvRadius",
     "The inverse hadron radius squared used in the overlap function",
     &MPIHandler::invRadius_, GeV2, 2.0*GeV2, 0.2*GeV2, 4.0*GeV2,
     true, false, Interface::limited);

  static ParVector<MPIHandler,int> interfaceadditionalMultiplicities
    ("additionalMultiplicities",
     "specify the multiplicities of secondary hard processes (multiple parton scattering)",
     &MPIHandler::additionalMultiplicities_, 
     -1, 0, 0, 3,
     false, false, true);

  static Parameter<MPIHandler,int> interfaceIdenticalToUE
    ("IdenticalToUE",
     "Specify which of the hard processes is identical to the UE one (QCD dijets)",
     &MPIHandler::identicalToUE_, -1, 0, 0,
     false, false, Interface::nolimits);

  static Parameter<MPIHandler,Energy> interfacePtOfQCDProc
    ("PtOfQCDProc",
     "Specify the value of the pt cutoff for the process that is identical to the UE one",
     &MPIHandler::PtOfQCDProc_, GeV, -1.0*GeV, ZERO, ZERO,
     false, false, Interface::nolimits);

  static Parameter<MPIHandler,double> interfacecolourDisrupt
    ("colourDisrupt",
     "Fraction of connections to additional subprocesses, which are colour disrupted.",
     &MPIHandler::colourDisrupt_, 
     0.0, 0.0, 1.0, 
     false, false, Interface::limited);

  
  static Switch<MPIHandler,bool> interfacesoftInt
    ("softInt",
     "Switch to enable soft interactions",
     &MPIHandler::softInt_, true, false, false);

  static SwitchOption interfacesoftIntYes
    (interfacesoftInt,
     "Yes",
     "enable the two component model",
     true);
  static SwitchOption interfacesoftIntNo
    (interfacesoftInt,
     "No",
     "disable the model",
     false);


  static Switch<MPIHandler,unsigned int> interEnergyExtrapolation
    ("EnergyExtrapolation",
     "Switch to ignore the cuts object at MPIHandler:Cuts[0]. "
     "Instead, extrapolate the pt cut.",
     &MPIHandler::energyExtrapolation_, 2, false, false);
  static SwitchOption interEnergyExtrapolationLog
    (interEnergyExtrapolation,
     "Log",
     "Use logarithmic dependence, ptmin = A * log (sqrt(s) / B).",
     1);
  static SwitchOption interEnergyExtrapolationPower
    (interEnergyExtrapolation,
     "Power",
     "Use power law, ptmin = pt_0 * (sqrt(s) / E_0)^b.",
     2);
  static SwitchOption interEnergyExtrapolationNo
    (interEnergyExtrapolation,
     "No",
     "Use manually set value for the minimal pt, "
     "specified in MPIHandler:Cuts[0]:OneCuts[0]:MinKT.",
     0);

  static Parameter<MPIHandler,Energy> interfaceEEparamA
    ("EEparamA",
     "Parameter A in the empirical parametrization "
     "ptmin = A * log (sqrt(s) / B)",
     &MPIHandler::EEparamA_, GeV, 0.6*GeV, ZERO, Constants::MaxEnergy,
     false, false, Interface::limited);

  static Parameter<MPIHandler,Energy> interfaceEEparamB
    ("EEparamB",
     "Parameter B in the empirical parametrization "
     "ptmin = A * log (sqrt(s) / B)",
     &MPIHandler::EEparamB_, GeV, 39.0*GeV, ZERO, Constants::MaxEnergy,
     false, false, Interface::limited);

  static Switch<MPIHandler,bool> interfacetwoComp
    ("twoComp",
     "switch to enable the model with a different radius for soft interactions",
     &MPIHandler::twoComp_, true, false, false);

  static SwitchOption interfacetwoCompYes
    (interfacetwoComp,
     "Yes",
     "enable the two component model",
     true);
  static SwitchOption interfacetwoCompNo
    (interfacetwoComp,
     "No",
     "disable the model",
     false);


  static Parameter<MPIHandler,CrossSection> interfaceMeasuredTotalXSec
    ("MeasuredTotalXSec",
     "Value for the total cross section (assuming rho=0). If non-zero, this "
     "overwrites the Donnachie-Landshoff parametrizations.",
     &MPIHandler::totalXSecExp_, millibarn, 0.0*millibarn, 0.0*millibarn, 0*millibarn,
     false, false, Interface::lowerlim);

  
  static Switch<MPIHandler,unsigned int> interfaceDLmode
    ("DLmode",
     "Choice of Donnachie-Landshoff parametrization for the total cross section.",
     &MPIHandler::DLmode_, 2, false, false);
  static SwitchOption interfaceDLmodeStandard
    (interfaceDLmode,
     "Standard",
     "Standard parametrization with s**0.08",
     1);
  static SwitchOption interfaceDLmodeCDF
    (interfaceDLmode,
     "CDF",
     "Standard parametrization but normalization fixed to CDF's measured value",
     2);
  static SwitchOption interfaceDLmodeNew
    (interfaceDLmode,
     "New",
     "Parametrization taking hard and soft pomeron contributions into account",
     3);

  static Parameter<MPIHandler,Energy> interfaceReferenceScale
    ("ReferenceScale",
     "The reference energy for power law energy extrapolation of pTmin",
     &MPIHandler::refScale_, GeV, 7000.0*GeV, 0.0*GeV, 20000.*GeV,
     false, false, Interface::limited);

  static Parameter<MPIHandler,Energy> interfacepTmin0
    ("pTmin0",
     "The pTmin at the reference scale for power law extrapolation of pTmin.",
     &MPIHandler::pT0_, GeV, 3.11*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<MPIHandler,double> interfacePower
    ("Power",
     "The power for power law extrapolation of the pTmin cut-off.",
     &MPIHandler::b_, 0.21, 0.0, 10.0,
     false, false, Interface::limited);
}
