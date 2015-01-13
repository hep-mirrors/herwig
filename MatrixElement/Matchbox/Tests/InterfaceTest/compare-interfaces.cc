  // This program uses many HepMC events and compares them.
  //
  // @author Johannes Bellm <johannes.bellm@kit.de>
  //

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/CompareGenEvent.h"
#include <gsl/gsl_randist.h>
#include <vector>
#include <iostream>     // std::cout, std::endl
#include <iomanip>

using namespace std;
using namespace HepMC;

  // check if a file exists

bool fileexists(string file) {
  bool result = false;
  fstream fin;
  fin.open(file.c_str(), ios::in);
  if (fin.is_open()) {
    result = true;
  }
  fin.close();
  return result;
}

void run() {
  
  
  
    // Specify the input files and the ascii output file
  std::vector<string> generators;
  if (fileexists("Truth.hepmc"))generators.push_back("Truth");
  if (fileexists("Buildin.hepmc"))generators.push_back("Buildin");
  if (fileexists("MadGraph.hepmc"))generators.push_back("MadGraph");
  if (fileexists("OpenLoops.hepmc"))generators.push_back("OpenLoops");
  if (fileexists("GoSam.hepmc"))generators.push_back("GoSam");
  if (fileexists("njet.hepmc"))generators.push_back("njet");
  
  
  if (generators.size()==0) {
    cout<< "No hepmc files found.";
    return;
  }
  if (generators.size()==1) {
    cout<< "Only one generator found. No comparison possible.";
    return;
  }
  cout<<"\nFound "<<generators.size()<<" Generators"<<flush;
  
  std::vector<HepMC::IO_GenEvent*> EventStreams;
  
  for (std::vector<string>::iterator it = generators.begin();it!=generators.end();it++){
    HepMC::IO_GenEvent* tmp= new HepMC::IO_GenEvent((*it+".hepmc"), ios::in);
    EventStreams.push_back(tmp);
  }
  
  
  std::vector<HepMC::GenEvent*> Events;
  bool stillevents=true;
    // get the first signal event
  for (std::vector<HepMC::IO_GenEvent*>::iterator it = EventStreams.begin();it!=EventStreams.end();it++){
    HepMC::GenEvent* evt =  (*it)->read_next_event();
    if (!evt)stillevents&=false;
    Events.push_back(evt);
  }
  if (!stillevents) {
    cout<<"An error accured  with the event files";
  }
  
  
    // loop until we run out of events
  while (stillevents) {
    std::vector<string>::iterator gen= generators.begin();
    std::vector<HepMC::GenEvent*>::iterator ev= Events.begin();
    cout<<"\n"<<*gen<<setw(10-(*gen).size())<<" ";
    cout<< std::setprecision(9) <<setw(14)<<((**(Events.begin())).weights()[0]);
    gen++;
    ev++;
    for(;ev!= Events.end();ev++,gen++){
      cout <<"  "<< *gen<<"/"<<generators[0]<<"=" <<setw(8)<< std::setprecision(6)<<((**ev).weights()[0])/((**(Events.begin())).weights()[0])<<flush;
      delete (*ev);
    }
    
      // read the next events
    std::vector<HepMC::GenEvent*>::iterator it2=Events.begin();
    for (std::vector<HepMC::IO_GenEvent*>::iterator it = EventStreams.begin();it!=EventStreams.end();it++,it2++){
      
      (**it)>>(*it2);
      if (!(*it2))stillevents&=false;
    }
  }
  
  cout<<"\nEnd of comparison.\n\n ";
}

int main() {
  run();
  return 0;
}

