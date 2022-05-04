// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ExperimentalOmnesFunction class.
//

#include "ExperimentalOmnesFunction.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ExperimentalOmnesFunction::ExperimentalOmnesFunction() {
  // initialization of the experimental function
  initialize_ =false;
  nPoints_=100;
  energy_ = {300*MeV, 320*MeV, 340*MeV, 360*MeV, 380*MeV,
	     400*MeV, 420*MeV, 440*MeV, 460*MeV, 480*MeV,
	     500*MeV, 520*MeV, 540*MeV, 560*MeV, 580*MeV,
	     600*MeV, 620*MeV, 640*MeV, 660*MeV, 680*MeV,
	     700*MeV, 720*MeV, 740*MeV, 760*MeV, 780*MeV,
	     800*MeV, 820*MeV, 840*MeV, 860*MeV, 880*MeV,
	     900*MeV, 920*MeV, 940*MeV, 960*MeV, 980*MeV};
  phase_={  00.1,     0.4,     0.7,     1.0,     1.5,
	    02.0,     2.5,     3.2,     4.0,     4.9,
	    05.9,     7.1,     8.5,    10.1,    12.1,
	    14.4,    17.3,    20.9,    25.4,    31.2,
	    38.7,    48.4,    60.6,    74.9,    90.0, 
	    103.8,   115.3,   124.3,   131.3,   136.7,
	    141.0,   144.5,   147.3,   149.7,   151.8};
  // experimental omnes function 
  omnesEnergy_ ={ 282.534*MeV, 289.32 *MeV, 296.106*MeV, 302.893*MeV, 309.679*MeV,
		  316.466*MeV, 323.252*MeV, 330.038*MeV, 336.825*MeV, 343.611*MeV,
		  350.398*MeV, 357.184*MeV, 363.97 *MeV, 370.757*MeV, 377.543*MeV,
		  384.33 *MeV, 391.116*MeV, 397.902*MeV, 404.689*MeV, 411.475*MeV,
		  418.261*MeV, 425.048*MeV, 431.834*MeV, 438.621*MeV, 445.407*MeV, 
		  452.193*MeV, 458.98 *MeV, 465.766*MeV, 472.553*MeV, 479.339*MeV,
		  486.125*MeV, 492.912*MeV, 499.698*MeV, 506.485*MeV, 513.271*MeV,
		  520.057*MeV, 526.844*MeV, 533.63 *MeV, 540.417*MeV, 547.203*MeV,
		  553.989*MeV, 560.776*MeV, 567.562*MeV, 574.349*MeV, 581.135*MeV,
		  587.921*MeV, 594.708*MeV, 601.494*MeV, 608.281*MeV, 615.067*MeV,
		  621.853*MeV, 628.64 *MeV, 635.426*MeV, 642.213*MeV, 648.999*MeV,
		  655.785*MeV, 662.572*MeV, 669.358*MeV, 676.145*MeV, 682.931*MeV,
		  689.717*MeV, 696.504*MeV, 703.29 *MeV, 710.077*MeV, 716.863*MeV,
		  723.649*MeV, 730.436*MeV, 737.222*MeV, 744.009*MeV, 750.795*MeV,
		  757.581*MeV, 764.368*MeV, 771.154*MeV, 777.94 *MeV, 784.727*MeV, 
		  791.513*MeV, 798.3  *MeV, 805.086*MeV, 811.872*MeV, 818.659*MeV, 
		  825.445*MeV, 832.232*MeV, 839.018*MeV, 845.804*MeV, 852.591*MeV, 
		  859.377*MeV, 866.164*MeV, 872.95 *MeV, 879.736*MeV, 886.523*MeV,
		  893.309*MeV, 900.096*MeV, 906.882*MeV, 913.668*MeV, 920.455*MeV, 
		  927.241*MeV, 934.028*MeV, 940.814*MeV, 947.6  *MeV, 954.387*MeV};
  omnesFunctionRe_ = { 0.860676  , 0.851786  , 0.843688  , 0.835827  , 0.828031  ,
		       0.820229  , 0.812370  , 0.804424  , 0.796354  , 0.788143  ,
		       0.779698  , 0.770939  , 0.761692  , 0.752707  , 0.743823  ,
		       0.735004  , 0.726091  , 0.717047  , 0.707862  , 0.698439  ,
		       0.688685  , 0.678510  , 0.668518  , 0.658481  , 0.648344  ,
		       0.638219  , 0.627989  , 0.617603  , 0.607222  , 0.596711  ,
		       0.586026  , 0.575280  , 0.564282  , 0.553067  , 0.541923  ,
		       0.530574  , 0.519112  , 0.507690  , 0.496055  , 0.484313  ,
		       0.472496  , 0.460245  , 0.447943  , 0.435766  , 0.423390  ,
		       0.410997  , 0.398510  , 0.385479  , 0.372458  , 0.359520  ,
		       0.346129  , 0.332837  , 0.319623  , 0.305858  , 0.292238  ,
		       0.278690  , 0.264391  , 0.250316  , 0.236400  , 0.221655  ,
		       0.207196  , 0.192956  , 0.177745  , 0.162833  , 0.148209  ,
		       0.132603  , 0.117202  , 0.102090  , 0.0862283 , 0.0703392 ,     
		       0.0545317 , 0.0383762 , 0.0219486 , 0.00518648,-0.0113217 ,    
		       -0.0280201 ,-0.045445  ,-0.0625479 ,-0.079748  ,-0.0978819 ,    
		       -0.11569   ,-0.133447  ,-0.152117  ,-0.170608  ,-0.189137  ,     
		       -0.208597  ,-0.227864  ,-0.247185  ,-0.267306  ,-0.287382  ,
		       -0.307707  ,-0.328882  ,-0.350103  ,-0.37178   ,-0.394464  ,
		       -0.417228  ,-0.440561  ,-0.464976  ,-0.490278  ,-0.517527};
  omnesFunctionIm_ = { 0.00243346, 0.000894972,-0.000612496,-0.00209178,-0.00354344,
		       -0.00496737,-0.00636316 ,-0.00773022 ,-0.00906769,-0.0103569 ,
		       -0.0116108 ,-0.0128658  ,-0.0145424  ,-0.0165746 ,-0.0186438 ,
		       -0.0206363 ,-0.0225379  ,-0.0243827  ,-0.0261488 ,-0.0278572 ,
		       -0.0295317 ,-0.0316349  ,-0.0339321  ,-0.0362345 ,-0.0386555 ,
		       -0.0410799 ,-0.0434534  ,-0.0459509  ,-0.0484302 ,-0.0508376 ,
		       -0.0533398 ,-0.0557937  ,-0.0581587  ,-0.0608612 ,-0.0635382 ,
		       -0.0661231 ,-0.068983   ,-0.0717604  ,-0.0744215 ,-0.0772635 ,
		       -0.0799845 ,-0.0825991  ,-0.0857537  ,-0.0888139 ,-0.0917441 ,
		       -0.0948263 ,-0.0977055  ,-0.100462   ,-0.103773  ,-0.106912  ,
		       -0.109931  ,-0.113413   ,-0.116647   ,-0.119722  ,-0.123282  ,
		       -0.126521  ,-0.129593   ,-0.133324   ,-0.136691  ,-0.139854  ,
		       -0.143729  ,-0.14718    ,-0.150356   ,-0.154353  ,-0.157926  ,
		       -0.161133  ,-0.165174   ,-0.168899   ,-0.172212  ,-0.176116  ,
		       -0.179892  ,-0.183445   ,-0.187134   ,-0.190947  ,-0.195144  ,
		       -0.198771  ,-0.202443   ,-0.206906   ,-0.210561  ,-0.214207  ,
		       -0.218943  ,-0.222806   ,-0.226551   ,-0.231273  ,-0.235267  ,
		       -0.239178  ,-0.244082   ,-0.24836    ,-0.252492  ,-0.257394  ,
		       -0.261812  ,-0.266156   ,-0.271161   ,-0.275849  ,-0.280675  ,
		       -0.286275  ,-0.291716   ,-0.297353   ,-0.303621  ,-0.310452  };
  // integration cut parameter
  epsCut_=0.4*MeV;
  // size of the arrays
  nsizea_ = energy_.size(); nsizeb_ = omnesEnergy_.size();
}


IBPtr ExperimentalOmnesFunction::clone() const {
  return new_ptr(*this);
}

IBPtr ExperimentalOmnesFunction::fullclone() const {
  return new_ptr(*this);
}

void ExperimentalOmnesFunction::doinit() {
  using Constants::pi;
  OmnesFunction::doinit();
  if(initialize_) {
    // convert the phase shift into radians
    vector<double> radphase;
    for(unsigned int ix=0;ix<phase_.size();++ix) {
      radphase.push_back(phase_[ix]/180.*Constants::pi);
    }
    // set up an interpolator for this
    interpolator_ = make_InterpolatorPtr(radphase,energy_,3);
    // limits and step sizes
    Energy mPi = getParticleData(ParticleID::piplus)->mass();
    Energy mEta(getParticleData(ParticleID::etaprime)->mass());
    Energy moff(2.*mPi),upp(mEta),step((mEta-moff)/nPoints_);
    // intergrator
    GaussianIntegrator integrator;
    // integrand
    // loop for integrals
    double D1real,D1imag;
    Complex ii(0.,1.),answer;
    moff+=0.5*step;
    omnesFunctionRe_.clear();
    omnesFunctionIm_.clear();
    omnesEnergy_.clear();
    for( ;moff<upp;moff+=step) {
      s_ = sqr(moff);
      // piece between 0 and 1 GeV
      Energy2 moff2(sqr(moff)),eps2(sqr(epsCut_));
      D1real=-moff2*(integrator.value(*this,4.*sqr(mPi),moff2-eps2)+
		     integrator.value(*this,moff2+eps2,upp*upp))/pi;
      D1imag=-(*interpolator_)(moff);
      // piece above 1 GeV
      D1real+=-(*interpolator_)(upp)/pi*log(upp*upp/(upp*upp-moff*moff));
      // calculate the answer
      answer = exp(D1real+ii*D1imag);
      // put into the arrays
      omnesFunctionRe_.push_back(answer.real());
      omnesFunctionIm_.push_back(answer.imag());
      omnesEnergy_.push_back(moff);
    }
  }
}

void ExperimentalOmnesFunction::persistentOutput(PersistentOStream & os) const {
  os << ounit(energy_,MeV) << ounit(omnesEnergy_,MeV) 
     << phase_ << omnesFunctionRe_ << omnesFunctionIm_ << initialize_
     << nPoints_ << ounit(epsCut_,MeV);
}

void ExperimentalOmnesFunction::persistentInput(PersistentIStream & is, int) {
  is >> iunit(energy_,MeV) >> iunit(omnesEnergy_,MeV) 
     >> phase_ >> omnesFunctionRe_ >> omnesFunctionIm_ >> initialize_
     >> nPoints_ >> iunit(epsCut_,MeV);
}

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<ExperimentalOmnesFunction,OmnesFunction>
describeHerwigExperimentalOmnesFunction("Herwig::ExperimentalOmnesFunction", "HwFormFactors.so");
HERWIG_INTERPOLATOR_CLASSDESC(ExperimentalOmnesFunction,double,Energy)

void ExperimentalOmnesFunction::Init() {

  static ClassDocumentation<ExperimentalOmnesFunction> documentation
    ("The ExperimentalOmnesFunction class implements the Omnes function"
     " integrating the experimental measurement of the phase");

  static Switch<ExperimentalOmnesFunction,bool> interfaceInitializeOmnes
    ("Initialize",
     "Initialize the experimental version of the Omnes function.",
     &ExperimentalOmnesFunction::initialize_, false, false, false);
  static SwitchOption interfaceInitializeOmnesInitialize
    (interfaceInitializeOmnes,
     "Yes",
     "Perform the initialization",
     true);
  static SwitchOption interfaceInitializeOmnesNoInitialization
    (interfaceInitializeOmnes,
     "No",
     "No initialization",
     false);

  static ParVector<ExperimentalOmnesFunction,Energy> interfacePhase_Energy
    ("Phase_Energy",
     "The energy values for the phase shift for the experimental version of the"
     " Omnes function",
     &ExperimentalOmnesFunction::energy_, MeV, -1, 1.0*MeV, 300.0*MeV, 2000.0*MeV,
     false, false, true);

  static ParVector<ExperimentalOmnesFunction,double> interfacePhase_Shift
    ("Phase_Shift",
     "The experimental values of the phase shift for the experimental version"
     " of the Omnes function",
     &ExperimentalOmnesFunction::phase_, 1.0, -1, 0.0, 0.0, 1000.0,
     false, false, true);

  static ParVector<ExperimentalOmnesFunction,Energy> interfaceOmnesEnergy
    ("OmnesEnergy",
     "The energy values for the interpolation of the experimental Omnes function",
     &ExperimentalOmnesFunction::omnesEnergy_, MeV, -1, 1.*MeV, 250.0*MeV, 2000.*MeV,
     false, false, true);

  static ParVector<ExperimentalOmnesFunction,double> interfaceOmnesReal
    ("OmnesReal",
     "The real part of the experimental Omnes function for the interpolation.",
     &ExperimentalOmnesFunction::omnesFunctionRe_, 1., -1, 1., -100.,
     100.,
     false, false, true);

  static ParVector<ExperimentalOmnesFunction,double> interfaceOmnesImag
    ("OmnesImag",
     "The imaginary part of the experimental Omnes function for the interpolation.",
     &ExperimentalOmnesFunction::omnesFunctionIm_, 1., -1, 1., -100.,
     100.,
     false, false, true);

  static Parameter<ExperimentalOmnesFunction,unsigned int> interfaceOmnesPoints
    ("OmnesPoints",
     "The number of points for the interpolation table for the experimental"
     " Omnes function.",
     &ExperimentalOmnesFunction::nPoints_, 100, 50, 200,
     false, false, true);

  static Parameter<ExperimentalOmnesFunction,Energy> interfaceOmnesCut
    ("OmnesCut",
     "The cut parameter for the integral in the experimental Omnes function.",
     &ExperimentalOmnesFunction::epsCut_, MeV, 0.1*MeV, 0.001*MeV, 1.0*MeV,
     false, false, true);

}

Complex ExperimentalOmnesFunction::D(Energy2 s) const {
  if(!oRe_) {
    oRe_ = make_InterpolatorPtr(omnesFunctionRe_,omnesEnergy_,3);
    oIm_ = make_InterpolatorPtr(omnesFunctionIm_,omnesEnergy_,3);
  }
  Energy q(sqrt(s)); Complex ii(0.,1.);
  return (*oRe_)(q)+ii*(*oIm_)(q);
}

void ExperimentalOmnesFunction::dataBaseOutput(ofstream & output,bool header,
					   bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::ExperimentalOmnesFunction " << name() << " \n";
  output << "newdef " << name() << ":Initialize " << initialize_      << "\n";
  output << "newdef " << name() << ":OmnesPoints     " << nPoints_         << "\n";
  output << "newdef " << name() << ":OmnesCut        " << epsCut_/MeV      << "\n";
  for(unsigned int ix=0;ix<energy_.size();++ix) {
    if(ix<nsizea_) {
      output << "newdef " << name() << ":Phase_Energy " << ix << "  " 
  	     << energy_[ix]/MeV << "\n";
      output << "newdef " << name() << ":Phase_Shift  " << ix << "  " 
  	     << phase_[ix]  << "\n";
    }
    else {
      output << "insert " << name() << ":Phase_Energy " << ix << "  " 
  	     << energy_[ix]/MeV << "\n";
      output << "insert " << name() << ":Phase_Shift  " << ix << "  " 
  	     << phase_[ix]  << "\n";
    }
  }
  for(unsigned int ix=0;ix<omnesEnergy_.size();++ix) {
    if(ix<nsizeb_) {
      output << "newdef " << name() << ":OmnesEnergy " << ix << "  " 
	     << omnesEnergy_[ix]/MeV << "\n";
      output << "newdef " << name() << ":OmnesReal " << ix << "  " 
	     << omnesFunctionRe_[ix] << "\n";
      output << "newdef " << name() << ":OmnesImag " << ix << "  " 
	     << omnesFunctionIm_[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":OmnesEnergy " << ix << "  " 
	     << omnesEnergy_[ix]/MeV << "\n";
      output << "insert " << name() << ":OmnesReal " << ix << "  " 
	     << omnesFunctionRe_[ix] << "\n";
      output << "insert " << name() << ":OmnesImag " << ix << "  " 
	     << omnesFunctionIm_[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
