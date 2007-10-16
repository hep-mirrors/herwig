// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MPIHandler class.
//

#include "MPIHandler.h"
#include "MPISampler.h"

#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Handlers/SubProcessHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/Handlers/LuminosityFunction.h"
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

/** Typedef for the MPISampler class.*/
typedef Ptr< MPISampler >::transient_pointer tMPISamplerPtr;

MPIHandler::MPIHandler()
  : theBinStrategy(2), theJmueo(1) {}

MPIHandler::MPIHandler(const MPIHandler & x)
  : Interfaced(x), LastXCombInfo<>(x), 
    theSampler(x.theSampler), theHandler(x.theHandler), 
    theCuts(x.theCuts), theSubProcesses(x.theSubProcesses),
    theXCombs(x.theXCombs), theXSecs(x.theXSecs),
    theBinStrategy(x.theBinStrategy), theMEXMap(x.theMEXMap),
    theMaxDims(x.theMaxDims), theMultiplicities(x.theMultiplicities),
    theJmueo(x.theJmueo), theRadius(x.theRadius) {}


MPIHandler::~MPIHandler() {}

void MPIHandler::initialize() {
  
  Energy maxEnergy = lumiFn().maximumCMEnergy();

  xCombs().clear();
  xSecs().clear();

  cuts()->initialize(sqr(maxEnergy), lumiFn().Y());

  for ( SubHandlerList::const_iterator sit = subProcesses().begin();
        sit != subProcesses().end(); ++sit ) {

    CutsPtr kincuts = (**sit).cuts()? (**sit).cuts(): cuts();
    if ( (**sit).cuts() ) kincuts->initialize(sqr(maxEnergy), lumiFn().Y());
    PExtrPtr pextract = (**sit).pExtractor();

    // Use an empty ckkw handler for the additional interactions:
    tCascHdlPtr ckkw = tCascHdlPtr();

    PartonPairVec vpc = pextract->getPartons(maxEnergy, incoming(), *kincuts);

    // The last parton bin pair was in fact the bins corresponding to
    // the incoming particles, so we remove them, but save them to
    // keep them referenced.
    PBPair orig = vpc.back();
    vpc.pop_back();

    for ( PartonPairVec::iterator ppit = vpc.begin();
          ppit != vpc.end(); ++ppit )
      for ( MEVector::const_iterator meit = (**sit).MEs().begin();
            meit != (**sit).MEs().end(); ++meit )
	addME(maxEnergy, *sit, pextract, kincuts, ckkw, *meit, *ppit);
  }

  xSecs().resize(xCombs().size());

  theMaxDims.clear();
  switch ( binStrategy() ) {
  case 0: {
    theMaxDims.push_back(0);
    for ( int i = 0, N = xCombs().size(); i < N; ++i )
      theMaxDims[0] = max(theMaxDims[0], xCombs()[i]->nDim());
    break;
  }
  case 1: {
    for ( int i = 0, N = xCombs().size(); i < N; ++i )
      theMEXMap[xCombs()[i]->matrixElement()].push_back(xCombs()[i]);
    MEXMap::const_iterator mei = theMEXMap.begin();
    for ( int i = 0, N = theMEXMap.size(); i < N; ++i, ++mei) {
      theMaxDims.push_back(0);
      for ( int j = 0, M = mei->second.size(); j < M; ++j )
        theMaxDims[i] = max(theMaxDims[i], mei->second[j]->nDim());
    }
    break;
  }
  case 2: {
    for ( int i = 0, N = xCombs().size(); i < N; ++i )
      theMaxDims.push_back(xCombs()[i]->nDim());
    break;
  }
  }

  tMPISamplerPtr smplr = dynamic_ptr_cast<tMPISamplerPtr>( sampler() );
  smplr->setMPIHandler(this);
  smplr->initialize();

  for ( int i = 0, N = xCombs().size(); i < N; ++i )
    xCombs()[i]->reset();

  double weight(0);
  //sample N PSpoints to get an estimate of the xsec
  for(unsigned int i=0; i<1000; i++){
    weight = sampler()->generate();
    tStdXCombPtr lastXC = select(sampler()->lastBin(), weight);
    weight /= lastXC->matrixElement()->preWeight();
  }

  ofstream file;
  file.open("multi.test");
  file.close();

  file.open("UE.out");
  Stat tot;

  statistics(file, tot);
  
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
  
  //now calculate the indivudual Probabilities
  XSVector UEXSecs;
  UEXSecs.push_back(tot.xSec());
  //  UEXSecs.push_back(99*millibarn);
  Probs(UEXSecs);
  UEXSecs.clear();

}

void MPIHandler::
addME(Energy maxEnergy, tSubHdlPtr sub, tPExtrPtr extractor, tCutsPtr cuts,
      tCascHdlPtr ckkw, tMEPtr me, const PBPair & pBins) {

  typedef MEBase::DiagramVector DiagramVector;
  typedef map<string,DiagramVector> DiagramMap;
  cPDPair pin(pBins.first->parton(), pBins.second->parton());
  DiagramVector diag = me->diagrams();
  DiagramMap tdiag;
  DiagramMap tmdiag;
  for ( int i = 0, N = diag.size(); i < N; ++i ) {
    if ( diag[i]->partons()[0] == pin.first &&
         diag[i]->partons()[1] == pin.second )
      tdiag[diag[i]->getTag()].push_back(diag[i]);
    if ( diag[i]->partons()[0] == pin.second &&
         diag[i]->partons()[1] == pin.first )
      tmdiag[diag[i]->getTag()].push_back(diag[i]);
  }

  bool mirror = false;
  if ( ( mirror = tdiag.empty() ) ) tdiag = tmdiag;
  for ( DiagramMap::iterator dit = tdiag.begin(); dit != tdiag.end(); ++dit ) {

    //todo: hope that it is no problem that I take the EventHandler here and not the MPIHandler:
    StdXCombPtr xcomb =
      new_ptr(StandardXComb(maxEnergy, incoming(), eventHandler(), 
			    sub, extractor, ckkw, pBins, cuts, me, dit->second, mirror));

    if ( xcomb->checkInit() ) xCombs().push_back(xcomb);

    else generator()->logWarning( 
      InitError() << "The matrix element '"
      << xcomb->matrixElement()->name() << "' cannot generate the diagram '"
      << dit->first << "' when used together with the parton extractor '"
      << xcomb->pExtractor()->name()
      << "'. The corresponding diagram is switched off." << Exception::warning);
  }
}

tStdXCombPtr MPIHandler::select(int bin, double weight) {

  int i = upper_bound(xSecs().begin(), xSecs().end(), UseRandom::rnd()*xSecs().back())
    - xSecs().begin();
  tStdXCombPtr lastXC;
  switch ( binStrategy() ) {
  case 0:
    lastXC = xCombs()[i];
    break;
  case 1: {
    MEXMap::iterator mei = theMEXMap.begin();
    for ( int j = 0; j < bin; ++j) ++mei;
    lastXC = mei->second[i];
    break;
  }
  case 2:
    lastXC = xCombs()[bin];
    break;
  }
  // clean up the old XComb object before switching to a new one
  if ( theLastXComb && theLastXComb != lastXC ) theLastXComb->clean();
  theLastXComb = lastXC;

  lastXC->select(weight);
  lastXC->accept();
  lastXC->matrixElement()->setXComb(lastXC);
  return lastXC;
}


CrossSection MPIHandler::
dSigDR(const pair<double,double> ll, Energy2 maxS,
       int ibin, int nr, const double * r) {

  PPair inc = make_pair(incoming().first->produceParticle(),
                        incoming().second->produceParticle());
  SimplePhaseSpace::CMS(inc, maxS);

  XVector xv;
  switch ( binStrategy() ) {
  case 0:
    xv = xCombs();
    break;
  case 1: {
    MEXMap::iterator mei = theMEXMap.begin();
    for ( int i = 0; i < ibin; ++i) ++mei;
    xv = mei->second;
    break;
  }
  case 2:
    xv = XVector(1, xCombs()[ibin]);
    break;
  }

  xSecs().resize(xv.size());
  for ( int i = 0, N = xv.size(); i < N; ++i ) xv[i]->prepare(inc);
  CrossSection sum = 0.0*nanobarn;
  for ( int i = 0, N = xv.size(); i < N; ++i )
    xSecs()[i] = ( sum += xv[i]->dSigDR(ll, nr, r) );

  return sum;
}


CrossSection MPIHandler::dSigDR(const vector<double> & r) {
  double jac = 1.0;
  pair<double,double> ll = lumiFn().generateLL(&r[0], jac);
  Energy2 maxS = sqr(lumiFn().maximumCMEnergy())/exp(ll.first + ll.second);
  int bin = sampler()->lastBin();
  CrossSection x = jac*lumiFn().value(incoming(), ll.first, ll.second)
    *dSigDR(ll, maxS, bin, nDim(bin) - lumiDim(), &r[lumiDim()]);
  return x;
}

int MPIHandler::nBins() const {
  switch ( binStrategy() ) {
  case 0: return 1;
  case 1: return theMEXMap.size();
  case 2: return xCombs().size();
  }
  return -1;
}

void MPIHandler::doinitrun() {
  
  Interfaced::doinitrun();
  sampler()->initrun();
  theHandler = generator()->currentEventHandler(); 
  //stop if the EventHandler is not present:
  assert(theHandler);
}


void MPIHandler::statistics(ostream & os, Stat & tot) const {

  if ( statLevel() == 0 ) return;
  map<cPDPair, Stat> partonMap;
  map<MEPtr, Stat> meMap;
  map<PExtrPtr, Stat> extractMap;
  //  Stat tot;

  for ( int i = 0, N = xCombs().size(); i < N; ++i ) {
    const StandardXComb & x = *xCombs()[i];
    Stat s;
    s = Stat(x.stats().attempts(), x.stats().accepted(),
             x.stats().sumWeights(), sampler()->integratedXSec(),
             sampler()->sumWeights());
    partonMap[x.partons()] += s;
    meMap[x.matrixElement()] += s;
    extractMap[x.pExtractor()] += s;
    tot += s;
  }

  string line = "======================================="
    "=======================================\n";

  if ( tot.accepted <= 0 ) {
    os << line << "No events generated by event handler '" << name() << "'."
       << endl;
    return;
  }

  os << line << "Statistics for event handler \'" << name() << "\':\n"
     << "                                       "
     << "generated    number of    Cross-section\n"
     << "                                       "
     << "   events     attempts             (nb)\n";

  os << line << "Total:" << setw(42) << tot.accepted << setw(13)
     << tot.attempted << setw(17) << tot.xSec()/nanobarn << endl
     << line;

  if ( statLevel() == 1 ) return;

  os << "Per matrix element breakdown:\n";
  for ( map<MEPtr, Stat>::iterator i = meMap.begin();
        i != meMap.end(); ++i ) {
    string n = i->first->name();
    n.resize(37, ' ');
    os << n << setw(11) << i->second.accepted << setw(13)
       << i->second.attempted << setw(17) << i->second.xSec()/nanobarn << endl;
  }
  os << line;

  if ( statLevel() == 2 ) return;

  os << "Per parton extractor breakdown:\n";
  for ( map<PExtrPtr, Stat>::iterator i = extractMap.begin();
        i != extractMap.end(); ++i ) {
    string n = i->first->name();
    n.resize(37, ' ');
    os << n << setw(11) << i->second.accepted << setw(13)
       << i->second.attempted << setw(17) << i->second.xSec()/millibarn << endl;
  }
  os << line;

  os << "Per incoming partons breakdown:\n";
  for ( map<cPDPair, Stat>::iterator i = partonMap.begin();
        i != partonMap.end(); ++i ) {
    string n = i->first.first->PDGName() + " " + i->first.second->PDGName();
    n.resize(37, ' ');
    os << n << setw(11) << i->second.accepted << setw(13)
       << i->second.attempted << setw(17) << i->second.xSec()/millibarn << endl;
  }
  os << line;

  if ( statLevel() == 3 ) return;

  os << "Detailed breakdown:\n";
  double xsectot = sampler()->integratedXSec()/
    (sampler()->sumWeights()*nanobarn);
  for ( int i = 0, N = xCombs().size(); i < N; ++i ) {
    const StandardXComb & x = *xCombs()[i];
    os << "(" << x.pExtractor()->name() << ") "
       << x.partons().first->PDGName() << " "
       << x.partons().second->PDGName()

       << " (" << x.matrixElement()->name() << " "
       << x.lastDiagram()->getTag() << ") " << endl
       << setw(48) << x.stats().accepted() << setw(13) << x.stats().attempts()
       << setw(17) << x.stats().sumWeights()*xsectot << endl;
  }

  os << line;
}

void MPIHandler::Probs(XSVector UEXSecs) {
  GaussianIntegrator integrator;
  ofstream file;
  file.open("multi.test", ios::app);
  file << "Probabilities:\n";
  
  unsigned int i(1);
  double P(0.0), AvgN(0.0);
  Length bmax(500.0*sqrt(millibarn));

  //currently only one UE process is possible so check that.
  assert(UEXSecs.size() == 1);

  for ( XSVector::const_iterator it = UEXSecs.begin();
        it != UEXSecs.end(); ++it ) {
    i = 1;
    Eikonalization inelint(this, *it, -1);//get the inel xsec
    do{
      //      cout << "debug: add integrand i = " << i << endl;
      Eikonalization integrand(this, *it, i);
      
      if(i>10) bmax = 10.0*sqrt(millibarn);
      if(theJmueo){
	P = integrator.value(integrand, Length(), bmax)/(*it);
      }else{
	P = integrator.value(integrand, Length(), bmax) /
	integrator.value(inelint, Length(), bmax);
      }
      AvgN += P*(i-1);
      
      file << "i = " << i-1 << ", P = " << P << "\n";
      //store the probability
      theMultiplicities.insert(P, i-1);

      i++;
    } while ( (i < 100) && (i < 5 || P > 1.e-15) );

    file << "------------------------------------------------\n";
    file << "AvgN: " << AvgN << endl;    
    file.close();
  }
  
}

double MPIHandler::factorial (unsigned int n) {
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
    if(theHandler->theJmueo)
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


void MPIHandler::persistentOutput(PersistentOStream & os) const {
  os << theSubProcesses << theCuts << theLastXComb
     << theXCombs << ounit(theXSecs, nanobarn)
     << theBinStrategy << theMaxDims << theMEXMap
     << theMultiplicities << theSampler << theHandler
     << theJmueo << ounit(theRadius, GeV2);
}

void MPIHandler::persistentInput(PersistentIStream & is, int) {
  is >> theSubProcesses >> theCuts >> theLastXComb
     >> theXCombs >> iunit(theXSecs, nanobarn)
     >> theBinStrategy >> theMaxDims >> theMEXMap
     >> theMultiplicities >> theSampler >> theHandler
     >> theJmueo >> iunit(theRadius, GeV2);
}

ClassDescription<MPIHandler> MPIHandler::initMPIHandler;
// Definition of the static class description member.

void MPIHandler::Init() {

  static ClassDocumentation<MPIHandler> documentation
    ("There is soon documentation for the MPIHandler class");

  static Reference<MPIHandler,SamplerBase> interfaceSampler
    ("Sampler",
     "The phase space sampler responsible for generating phase space"
     "points according to the cross section given by this handler",
     &MPIHandler::theSampler, false, false, true, true);

  static RefVector<MPIHandler,SubProcessHandler> interfaceSubhandlers
    ("SubProcessHandlers",
     "The list of sub-process handlers used in this EventHandler. ",
     &MPIHandler::theSubProcesses, 0, false, false, true, false);

  static Reference<MPIHandler,Cuts> interfaceCuts
    ("Cuts",
     "Common kinematical cuts for this MultipleInteractionHandler. These cuts "
     "should not be ovveridden in individual sub-process handlers.",
     &MPIHandler::theCuts, false, false, true, false);


  static Parameter<MPIHandler,Energy2> interfaceRadius
    ("Radius",
     "The inverse hadron radius squared, used in the overlap function",
     &MPIHandler::theRadius, GeV2, 0.71*GeV2, 0.0*GeV2, 0*GeV2,
     true, false, Interface::lowerlim);


  static Switch<MPIHandler,int> interfaceJmueo
    ("Jmueo",
     "This option determines in which mode the UE algorithm runs. "
     "0 for UE under low pt jets, 1 for efficient generation "
     "of UE activity with a rare signal process.",
     &MPIHandler::theJmueo, 1, false, false);

  static SwitchOption interfaceJmueo0
    (interfaceJmueo,
     "lowpt",
     "Signal process has similar cross section than UE.",
     0);

  static SwitchOption interfaceJmueo1
    (interfaceJmueo,
     "rare",
     "Signal process has a much smaller cross section than UE.",
     1);

  static Switch<MPIHandler,int> interfaceBinStrategy
    ("BinStrategy",
     "The strategy to be used when sampling different ThePEG::XComb "
     "objects. An ThePEG::XComb objet represents a pair of incoming "
     "parton types as defined by a ThePEG::PartonExtractor and a "
     "matrix element.",
     &MPIHandler::theBinStrategy, 2, false, false);

  static SwitchOption interfaceBinStrategy0
    (interfaceBinStrategy,
     "AllAtOnce",
     "All bins are sampled together.",
     0);

  static SwitchOption interfaceBinStrategy1
    (interfaceBinStrategy,
     "PerME",
     "All bins which have the same matrix element object are sampled together.",
     1);

  static SwitchOption interfaceBinStrategy2
    (interfaceBinStrategy,
     "Individual",
     "All bins are sampled individually.",
     2);


}

InvArea MPIHandler::OverlapFunction(Length b) {
  InvLength mu = sqrt(theRadius)/hbarc;
  return (sqr(mu)/96/Constants::pi)*pow(mu*b, 3)*(gsl_sf_bessel_Kn(3, mu*b));
}
