// -*- C++ -*-
//
// MPIHandler.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/Handlers/CascadeHandler.h"
#include "ThePEG/Cuts/Cuts.h"

#include "Herwig++/Utilities/GaussianIntegrator.h"

#include "gsl/gsl_sf_bessel.h"


#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MPIHandler.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MPIHandler::MPIHandler()
  : theAlgorithm(2) {}

MPIHandler::MPIHandler(const MPIHandler & x)
  : Interfaced(x), 
    theHandler(x.theHandler), theSubProcesses(x.theSubProcesses),
    theCuts(x.theCuts), theProcessHandlers(x.theProcessHandlers),
    theMultiplicities(x.theMultiplicities),
    theAlgorithm(x.theAlgorithm), theInvRadius(x.theInvRadius) {}


MPIHandler::~MPIHandler() {}

void MPIHandler::finalize() {
  if( beamOK() )
    statistics("UE.out");
}
void MPIHandler::initialize() {
  
  theHandler = generator()->currentEventHandler(); 
  //stop if the EventHandler is not present:
  assert(theHandler);

  //check if MPI is wanted
  if( !beamOK() ){
    generator()->log() << "You have requested multiple parton-parton scattering,\n"
		       << "but the model is not forseen for the setup you chose.\n" 
		       << "Events will be produced without MPI.\n";
    return;
  }


  if( subProcesses().size() != cuts().size() ) 
    throw Exception() << "MPIHandler::each SubProcess needs a Cuts Object"
		      << "ReferenceVectors are not equal in size"
		      << Exception::runerror;

  for(unsigned int i=0; i<cuts().size(); i++){
    theProcessHandlers.push_back(new_ptr(ProcessHandler()));
    processHandlers().back()->initialize(subProcesses()[i], 
					 cuts()[i], theHandler);
  }


  for(unsigned int i=0; i<cuts().size(); i++)
    processHandlers()[i]->initrun();

  /*
  if(Algorithm()==0){
  
    //check out the eikonalization -1=inelastic, -2=total xsec
    Eikonalization integrand(this, tot.xSec(), -1);
    Eikonalization integrand_tot(this, tot.xSec(), -2);
    GaussianIntegrator integrator;

    string line = "======================================="
      "=======================================\n";
  
    CrossSection inel(integrator.value(integrand, Length(), 1000.*sqrt(millibarn))), 
      total(integrator.value(integrand_tot, Length(), 1000.*sqrt(millibarn)));

    file << "\nEikonalization results:\n"
         << setw(79)
         << "Cross-section (mb)\n"
         << line << "Inelastic cross-section" << setw(55) 
         << inel/millibarn << endl
         << "Total pp->X cross-section" << setw(53)
         << total/millibarn << endl << line 
         << "Average number of MPI" << setw(57) << tot.xSec()/inel << endl;

    file.close();
  }
  
  */


  //now calculate the individual Probabilities
  XSVector UEXSecs;
  UEXSecs.push_back(processHandlers()[0]->integratedXSec());

  Probs(UEXSecs);
  UEXSecs.clear();
}


void MPIHandler::statistics(string os) const {
  ofstream file;
  file.open(os.c_str());

  for(unsigned int i=0; i<cuts().size(); i++){
    Stat tot;
    file << "Process " << i << ":\n";
    processHandlers()[i]->statistics(file, tot);
    file << "\n";
  }
  file.close();
}

void MPIHandler::Probs(XSVector UEXSecs) {
  GaussianIntegrator integrator;
  unsigned int i(1);
  double P(0.0), AvgN(0.0);
  Length bmax(500.0*sqrt(millibarn));

  //currently only one UE process is possible so check that.
  assert(UEXSecs.size() == 1);
  ofstream file;
  file.open("probs.test");
  file << "hard process xsec: "
       << dynamic_ptr_cast<tStdEHPtr>(eventHandler())->integratedXSec()/millibarn
       << endl;             

  file << "UE process[0] xsec: " 
       << UEXSecs.front()/millibarn << endl; 

  for ( XSVector::const_iterator it = UEXSecs.begin();
        it != UEXSecs.end(); ++it ) {
    i = 1;
    Eikonalization inelint(this, *it, -1);//get the inel xsec
    do{
      //      cout << "debug: add integrand i = " << i << endl;
      Eikonalization integrand(this, *it, i);
      
      if(i>10) bmax = 10.0*sqrt(millibarn);
      if(theAlgorithm > 0){
	P = integrator.value(integrand, Length(), bmax)/(*it);
      }else{
	P = integrator.value(integrand, Length(), bmax) /
	integrator.value(inelint, Length(), bmax);
      }
      AvgN += P*(i-1);
      //store the probability
      theMultiplicities.insert(P, i-1);
      file << i-1 << " " << P << endl;
      i++;
    } while ( (i < 100) && (i < 5 || P > 1.e-15) );
  }
  file.close();  
}


// calculate the integrand
Length Eikonalization::operator() (Length b) const {
  //fac is just: db^2=fac*db despite that large number
  unsigned int n(0);
  Length fac(Constants::twopi*b);
  CrossSection sigma(theUneikXSec);
  InvArea Ab(theHandler->OverlapFunction(b));

  //total cross section wanted
  if(theoption == -2) return 2 * fac * ( 1 - exp(-Ab*sigma / 2.) );

  //inelastic cross section
  if(theoption == -1) return   fac * ( 1 - exp(-Ab*sigma) );

  //P_n*sigma. Described in MPIHandler.h
  if(theoption > 0){
    n=theoption;
    if(theHandler->theAlgorithm > 0)
      return fac / theHandler->factorial(n-1) * pow(Ab*sigma, double(n)) 
	* exp(-Ab*sigma);
    else
      return fac / theHandler->factorial(n) * pow(Ab*sigma, double(n)) 
	* exp(-Ab*sigma);
  }else{
    throw Exception() << "Parameter theoption in Struct Eikonalization in " 
		      << "MPIHandler.cc has not allowed value"
                      << Exception::runerror;
    return 0.0*meter;
  }
}

double MPIHandler::factorial (unsigned int n) const {
  static unsigned int max(100);
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

  if(n > max) 
        throw Exception() << "MPIHandler::factorial called with too large argument"
                      << Exception::runerror;
  else
    return f[n];
}

InvArea MPIHandler::OverlapFunction(Length b) const {
  InvLength mu = sqrt(theInvRadius)/hbarc;
  return (sqr(mu)/96/Constants::pi)*pow(mu*b, 3)*(gsl_sf_bessel_Kn(3, mu*b));
}

double MPIHandler::poisson(Length b, CrossSection sigma, unsigned int N) const {
  return pow(OverlapFunction(b)*sigma, (double)N)/factorial(N)
    *exp(-OverlapFunction(b)*sigma);
}

void MPIHandler::persistentOutput(PersistentOStream & os) const {
  os << theMultiplicities << theHandler
     << theSubProcesses << theCuts << theProcessHandlers
     << theAlgorithm << ounit(theInvRadius, GeV2);
}

void MPIHandler::persistentInput(PersistentIStream & is, int) {
  is >> theMultiplicities >> theHandler
     >> theSubProcesses >> theCuts >> theProcessHandlers
     >> theAlgorithm >> iunit(theInvRadius, GeV2);
}

ClassDescription<MPIHandler> MPIHandler::initMPIHandler;
// Definition of the static class description member.

void MPIHandler::Init() {

  static ClassDocumentation<MPIHandler> documentation
    ("There is soon documentation for the MPIHandler class");

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
     &MPIHandler::theInvRadius, GeV2, 1.3*GeV2, 0.2*GeV2, 4.0*GeV2,
     true, false, Interface::limited);


  static Switch<MPIHandler,int> interfaceAlgorithm
    ("Algorithm",
     "This option determines in which mode the UE algorithm runs. "
     "0 for UE under low pt jets, 1 for UE(jets) under highpt jets, "
     "2 for efficient generation "
     "of UE activity with a rare signal process.",
     &MPIHandler::theAlgorithm, 2, false, false);

  static SwitchOption interfaceAlgorithm0
    (interfaceAlgorithm,
     "lowpt",
     "Signal process has similar cross section than UE.",
     0);

  static SwitchOption interfaceAlgorithm1
    (interfaceAlgorithm,
     "highpt",
     "Signal process has a much smaller cross section "
     "than UE, but the same ME's",
     1);

  static SwitchOption interfaceAlgorithm2
    (interfaceAlgorithm,
     "rare",
     "Signal process has a much smaller cross section "
     "than UE and is a different process.",
     2);
}
