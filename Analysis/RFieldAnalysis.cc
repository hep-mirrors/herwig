// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RFieldAnalysis class.
//

#include "RFieldAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "RFieldAnalysis.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

class Jet{

public:

  Jet(){
    ptsum = 0.0*GeV; 
    ptweighted_phi = 0.0*GeV; 
    ptweighted_eta = 0.0*GeV;
  }
  void add(const Lorentz5Momentum &p){
    ptweighted_phi += p.perp()*p.phi();
    ptweighted_eta += p.perp()*p.eta();
    ptsum += p.perp();
  } 
  bool isInCircle(const Lorentz5Momentum &p, double R){
    double deta = eta() - p.eta();
    double phi1 = phi();
    double phi2 = p.phi();
    double dphi = fabs( phi1 - phi2 );
    if (dphi > Constants::pi) dphi = fabs( dphi - 2*Constants::pi );
    if( (sqr(deta) + sqr(dphi)) < sqr(R) ) return true;
    return false;
  }
  Energy perp() const {return ptsum;}
  double phi() const {return ptweighted_phi/ptsum;}
  double eta() const {return ptweighted_eta/ptsum;}

private:
  Energy ptsum;
  Energy ptweighted_phi;
  Energy ptweighted_eta;
};

typedef multimap<double, tPPtr> ptlist;
typedef multimap<double, Jet> Jetptlist;

/**
 * Returns the list of jets 
 */
vector<Jet> getJets(vector<tPPtr> particles, double R=0.7){
  ptlist tracks;
  ptlist::iterator t1, t2;
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

  for(Jetptlist::iterator jit=tmp.begin(); jit!=tmp.end(); ++jit)
    result.push_back(jit->second);

//  cerr << "return result: " << result[0].perp()/GeV << endl;
  return result;
}

RFieldAnalysis::~RFieldAnalysis() {}

void RFieldAnalysis::analyze(tEventPtr event, long , int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  /** get the final-state particles */
  tPVector particles=event->getFinalState();
  tPVector selection;
  Lorentz5Momentum p;
  int nch(0);//number of all charged particles
  int nchTow(0), nchTrans(0), nchAway(0);
  double ptsumTow(0.0), ptsumTrans(0.0), ptsumAway(0.0);
  double dphi(0.0), pt1(0.0);
  
  for (tPVector::const_iterator pit = particles.begin(); pit != particles.end(); ++pit){
      /** Select only the charged particles  */
      if( ChargedSelector::Check(**pit) ){
          //CDF detector simulation:
          if(UseRandom::rnd() < 0.08) continue;

	  p = (**pit).momentum();
	  /**
	     Select the particles according the noted selection cuts and do the
	     analysis (towards away and transverse region) only on this set of particles.
	  */
	  if ( fabs( p.eta() ) < 1 && p.perp() > 0.5*GeV )
	      selection.push_back( *pit );
      
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
      vector<Jet>::const_iterator itr = jets.begin();

      pt1 = itr->perp()/GeV;//leading jet scalar ptsum of constituent particles
      
      //overflow or underflow
      if(pt1 < (double)thelow || pt1 > (double)theup) return;
      //get the bin number
      int id((int)(floor(pt1))-thelow);
      
      for (tPVector::const_iterator pit = selection.begin(); pit != selection.end(); ++pit){

	  p = (**pit).momentum();

	  dphi = fabs( p.phi() - itr->phi() )*180.0/Constants::pi;
	  if (dphi > 180) dphi = fabs( dphi - 360.0 );
      
          //Towards region
	  if( dphi < 60 ){
            nchTow++;
            ptsumTow += p.perp()/GeV;
          }
          //Transverse region
	  if( dphi > 60 && dphi < 120 ){
            nchTrans++;
            ptsumTrans += p.perp()/GeV;
	  }
          //Away region
	  if( dphi > 120 ){
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

