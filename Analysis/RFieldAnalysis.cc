// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RFieldAnalysis class.
//

#include "RFieldAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "Herwig++/Interfaces/KtJetInterface.h"
#include "KtJet/KtEvent.h"
#include "KtJet/KtLorentzVector.h"


#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "RFieldAnalysis.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

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
    to the reconstructed jet with larges pt and write all the particles that 
    are in there to the ROOT file in the same manner as above.
  */
  if(selection.size()){

      KtJetInterface jet;
      KtJet::KtEvent ev(jet.convert(selection), 4, 2, 2, 0.7);//pt recom scheme
      /** Get the leading jet (largest pt) */
      vector<KtJet::KtLorentzVector> jets = ev.getJetsPt();
      vector<KtJet::KtLorentzVector>::const_iterator itr = jets.begin();
  
      pt1 = itr->perp()/1000.0;//leading jet pt
      
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
  os << theShowerHandler << thelow << theup << theDir;
}

void RFieldAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> theShowerHandler >> thelow >> theup >> theDir;
}

ClassDescription<RFieldAnalysis> RFieldAnalysis::initRFieldAnalysis;
// Definition of the static class description member.

void RFieldAnalysis::Init() {

  static ClassDocumentation<RFieldAnalysis> documentation
    ("There is no documentation for the RFieldAnalysis class");

  
  static Reference<RFieldAnalysis,ShowerHandler> interfaceShowerHandler
    ("ShowerHandler",
     "A reference to the ShowerHandler",
     &RFieldAnalysis::theShowerHandler, true, false, true, false, false);


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

