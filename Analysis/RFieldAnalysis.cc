// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RFieldAnalysis class.
//

#include "RFieldAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

double mapPhiPI(const double dphi){//map phi to -pi,pi
  if (dphi > Constants::pi || dphi < -Constants::pi ){ 
    if(dphi > Constants::pi)
      return dphi - 2*Constants::pi;
    else
      return dphi + 2*Constants::pi;
  }else{
    return dphi;
  }
}

double phiDiffPI(const double phi1, const double phi2){
    double dphi = fabs( phi1 - phi2 );
    return fabs(mapPhiPI(dphi));
}

class Jet{

public:

typedef vector<pair<double, Energy> > phiptlist;

  Jet(){
    ptsum = 0.0*GeV; 
    ptweighted_eta = 0.0*GeV;
    calcphi = true;
  }
  void add(const Lorentz5Momentum &p){
    ptweighted_eta += p.perp()*p.eta();
    ptsum += p.perp();
    phipt.push_back(make_pair(p.phi(), p.perp()));
    calcphi = true;
  } 

  bool isInCircle(const Lorentz5Momentum &p, double R){
    double deta = eta() - p.eta();
    double phi1 = phi();
    double phi2 = p.phi();
    double dphi = phiDiffPI(phi1, phi2);
    /*
    cerr << "check dR for: eta: " << p.eta() << " pt: " << p.perp()/GeV
	 << " phi: " << p.phi() << endl;
    cerr << "check dR: " << sqrt(sqr(deta)+sqr(dphi)) << endl;
    */
    if( (sqr(deta) + sqr(dphi)) < sqr(R) ) return true;
    return false;
  }
  Energy perp() const {return ptsum;}
  double phi() {
    if(!calcphi) return _phi;

    double phi1( phipt.front().first );
    double phi = phi1;

    //cerr << "entering phi()....\n";
    
    //cerr << "seed: " << phi1 << " " << phipt.front().second/GeV << endl;

    phiptlist::const_iterator it = phipt.begin();
    it++;
    while(it != phipt.end()){
      phi += it->second/ptsum*mapPhiPI(it->first - phi1);
      //cerr << "loop: " << it->first << " " << it->second/GeV << endl;
      ++it;
    }
    //cerr << "result: " << phi << " " << mapPhiPI(phi) << endl;
    _phi = mapPhiPI(phi);
    calcphi = false;//avoid calculating phi twice
    return _phi;
  }
  double eta() const {return ptweighted_eta/ptsum;}

private:
  bool calcphi;
  double _phi;
  phiptlist phipt;
  Energy ptsum;
  Energy ptweighted_eta;
};

typedef multimap<double, tPPtr> ptlist;
typedef multimap<double, Jet> Jetptlist;

/**
 * Returns the list of jets 
 */
vector<Jet> getJets(vector<tPPtr> particles, double R=0.7){
  ptlist tracks;
  ptlist::iterator t1, t2, t3;
  Jetptlist tmp;
  vector<Jet> result;
  
  for(unsigned int i=0; i<particles.size(); i++)
    //sort in decreasing pt
    tracks.insert(make_pair(-particles[i]->momentum().perp()/GeV, particles[i]));
    
  /* output all particles to work with
  for(t1=tracks.begin(); t1!=tracks.end(); ++t1)
    cerr << t1->first << "--" << t1->second << endl;
  */

  while( !tracks.empty() ){
    t1 = tracks.begin();

    Jet curjet;
    curjet.add(t1->second->momentum());
    tracks.erase(t1);

    t2 = tracks.begin();
    while (t2 != tracks.end()) {
      if( curjet.isInCircle(t2->second->momentum(), R) ){
        curjet.add(t2->second->momentum());
	tracks.erase(t2);
      }
      ++t2;
    }    
    tmp.insert(make_pair(-curjet.perp()/GeV, curjet));
//    cerr << "Jet with ptsum: " << curjet.perp()/GeV << endl;
  }

  for(Jetptlist::iterator jit=tmp.begin(); jit!=tmp.end(); ++jit){
    result.push_back(jit->second);
    //    cerr << "Jet axis: pt: " << -jit->first << " phi: " << jit->second.phi() << endl;
  }
  return result;
}


/**
 * A somewhat more sophisticated fastCDFSim than using a finite
 * acceptance constant over pt. It is extracted from plot 19 in 
 * CDF/DOC/MIN BIAS/CDFR/8593
 * Track selection and counting efficiency in minimum bias, by
 * N. Moggi, M. Mussini, F. Rimondi
 */
bool keep(Energy ){
  double acceptance(1.0);
/* 
  if( pt < 0.5*GeV)
    acceptance /= 1.5;

  if( pt > 0.5*GeV && pt <= 1.0*GeV )
    acceptance /= 1.2;

  if( pt > 1.0*GeV && pt <= 1.5*GeV )
    acceptance /= 1.14;

  if( pt > 1.5*GeV && pt <= 2.0*GeV )
    acceptance /= 1.1;

  if( pt > 2.0*GeV && pt <= 2.5*GeV )
    acceptance /= 1.073;

  if( pt > 2.5*GeV && pt <= 3.0*GeV )
    acceptance /= 1.054;

  if( pt > 3.0*GeV && pt <= 4.0*GeV )
    acceptance /= 1.040;

  if( pt > 4.0*GeV && pt <= 5.0*GeV )
    acceptance /= 1.029;

  if( pt > 5.0*GeV )
    acceptance /= 1.02;
*/

  //use constant acceptance again
  acceptance = 0.92;

  if(UseRandom::rnd() < acceptance)
    return true;
  else
    return false;
}

RFieldAnalysis::RFieldAnalysis(): thelow(30), 
				  theup(50), theDir(".") {}


void RFieldAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  //equal size bins from thelow to theup
  for(int i=0; i<(theup-thelow); i++){
    theNTow.push_back(Statistic());
    theNTrans.push_back(Statistic());
    theNAway.push_back(Statistic());
    thePtsumTow.push_back(Statistic());
    thePtsumTrans.push_back(Statistic());
    thePtsumAway.push_back(Statistic());
  }
}

void RFieldAnalysis::analyze(tEventPtr event, long , int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  /** get the final-state particles */
  tPVector particles=event->getFinalState();
  tPVector selection;
  Lorentz5Momentum p;
  int nch(0);//number of all charged particles
  int nchTow(0), nchTrans(0), nchAway(0);
  double ptsumTow(0.0), ptsumTrans(0.0), ptsumAway(0.0);
  double dphi(0.0), pt1(0.0), philead(0.0);
  Energy pt(0.0*GeV);

  for (tPVector::const_iterator pit = particles.begin(); pit != particles.end(); ++pit){
      /** Select only the charged particles  */
      if( ChargedSelector::Check(**pit) ){

	  p = (**pit).momentum();
          pt = p.perp();
          //CDF detector simulation:
	  if(!keep(pt)) continue;
	  
	  /**
	     Select the particles according the noted selection cuts and do the
	     analysis (towards away and transverse region) only on this set of particles.
	  */
	  if ( fabs( p.eta() ) < 1 && pt > 0.5*GeV ){
	      selection.push_back( *pit );
	      //	      cerr << p.eta() << " " << pt/GeV << " " 
	      //		   << p.phi() << endl;
	  }      

	  nch++;//number of all charged particles
      }
  }
  /*
    get the "transverse" region, defined by the azimuthal angle difference
    to the reconstructed jet with largest pt.
  */
  if(selection.size()){

      vector<Jet> jets = getJets(selection);
      //Get the leading jet (largest pt)
      vector<Jet>::iterator itr = jets.begin();

      pt1 = itr->perp()/GeV;//leading jet scalar ptsum of constituent particles
      philead = itr->phi();

      //      cerr << "Hw-leadingjet: " << philead << " " << pt1 << endl;
      //overflow or underflow
      if(pt1 < (double)thelow || pt1 > (double)theup) return;
      //get the bin number
      int id((int)(floor(pt1))-thelow);

      //      cerr << "herwig: ljet " << pt1 << " " << itr->phi() << endl;
      //cerr << "-------------------------------------------\n";

      for (tPVector::const_iterator pit = selection.begin(); pit != selection.end(); ++pit){

	  p = (**pit).momentum();
	  dphi = phiDiffPI(p.phi(), philead);//*180.0/Constants::pi;
      
          //Towards region
	  if( dphi < Constants::pi/3.0 ){
            nchTow++;
            ptsumTow += p.perp()/GeV;
          }
          //Transverse region
	  if( dphi > Constants::pi/3.0 && dphi < 2.0*Constants::pi/3.0 ){
            nchTrans++;
            ptsumTrans += p.perp()/GeV;
	  }
          //Away region
	  if( dphi > 2.0*Constants::pi/3.0 ){
            nchAway++;
            ptsumAway += p.perp()/GeV;            
          }
      }
      theNTow[id]       += (double)nchTow;
      theNTrans[id]     += (double)nchTrans;
      theNAway[id]      += (double)nchAway;
      thePtsumTow[id]   += ptsumTow;
      thePtsumTrans[id] += ptsumTrans;
      thePtsumAway[id]  += ptsumAway;
  }
}

void RFieldAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  pair<double, int> tot = make_pair(0.0, 0);
  ofstream file;
  file.open("RField-observables.dat");
  file.close();

  generator()->log() << "\nChi^2 for R.Fields TVT analyses [Phys.Rev.D65:092002,2002]\n";

  chisq(theNTrans, theDir+"/nch_transv.dat", tot);
  chisq(theNTow, theDir+"/nch_tow.dat", tot);
  chisq(theNAway, theDir+"/nch_away.dat", tot);

  chisq(thePtsumTrans, theDir+"/ptsum_transv.dat", tot);
  chisq(thePtsumTow, theDir+"/ptsum_tow.dat", tot);
  chisq(thePtsumAway, theDir+"/ptsum_away.dat", tot); 

  generator()->log() << "Total Chi^2: " << tot.first << "/" << tot.second 
                     << " = " << tot.first/tot.second << endl;
}

void RFieldAnalysis::chisq(vector<Statistic> mc, string ascii, pair<double, int> &tot){

  vector<DataContainer *> data = getData(ascii);
  double datalow = data.front()->low;
  double datahigh= data.back()->up;
  ofstream file;
  file.open("RField-observables.dat", ios::app);
  //assume equal bin width:
  double width = (datahigh-datalow)/data.size();
  assert(width == 1.0);
  assert(datalow <= thelow && datahigh >= theup);
  
  int used_bins = (int)((theup-thelow)/width);
  int data_offset = (int)((thelow-datalow)/width);

  double chi2(0.0);
  double mean;

  file << "# histogram "+ascii
       << "\n# binlow binup mcmean mcerror\n";

  for(int i=0; i<used_bins; i++){
    mean = mc[i].mean(); 

    if(mean != 0)
      chi2 += sqr(mean - data[i+data_offset]->content) /
        (mc[i].mean_var() + sqr(data[i+data_offset]->error));

    file << data[i+data_offset]->low << " " << data[i+data_offset]->up
         << " " << mean << " " << mc[i].mean_stdDev() << endl;
  }
  file.close();
  generator()->log() << "Chi^2 for " << ascii << ":\n\t" << chi2 
                     << "/" << used_bins << " = " << chi2/used_bins << endl;
  
  tot.first += chi2;
  tot.second += used_bins;
  data.clear();
}

vector<DataContainer *> RFieldAnalysis::getData(string ascii){
  ifstream file(ascii.c_str(), ios::in);
  if(!file)
    throw Exception() << "Data file for Chi2 calculation " 
                      << "could not be opened in RFieldAnalysis.cc"
                      << Exception::runerror;

  DataContainer *data(0);
  char text[100] = "";
  vector<DataContainer *> table;

  //read the ascii file
  while( !file.eof() ){
    file.getline(text, 100);
    data = new DataContainer();
    sscanf(text,"%lf%lf%lf%lf", &data->low, &data->up, &data->content, &data->error);
    if((data->low || data->up) && data->content)
      table.push_back(data);
  }
  return table;
}

LorentzRotation RFieldAnalysis::transform(tEventPtr ) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void RFieldAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void RFieldAnalysis::analyze(tPPtr) {}

void RFieldAnalysis::persistentOutput(PersistentOStream & os) const {
  os << thelow << theup << theDir;
}

void RFieldAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> thelow >> theup >> theDir;
}

ClassDescription<RFieldAnalysis> RFieldAnalysis::initRFieldAnalysis;
// Definition of the static class description member.

void RFieldAnalysis::Init() {

  static ClassDocumentation<RFieldAnalysis> documentation
    ("There is no documentation for the RFieldAnalysis class");

  
  static Parameter<RFieldAnalysis,int> interfaceLowerBorder
    ("LowerBorder",
     "Lower border of chi^2 comparison to data.",
     &RFieldAnalysis::thelow, 30, 0, 0,
     true, false, Interface::lowerlim);

  static Parameter<RFieldAnalysis,int> interfaceUpperBorder
    ("UpperBorder",
     "Upper border of chi^2 comparison to data.",
     &RFieldAnalysis::theup, 50, 0, 0,
     true, false, Interface::lowerlim);

  
  static Parameter<RFieldAnalysis,string> interfaceDir
    ("Dir",
     "Directory where data files are stored",
     &RFieldAnalysis::theDir, ".",
     true, false);
}

