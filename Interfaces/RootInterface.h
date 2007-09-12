#ifndef ROOT_INTERFACE_H_
#define ROOT_INTERFACE_H_

// This is the declaration of the RootInterface class.
#include "Herwig++/Utilities/HerwigRun.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

namespace Herwig {

using namespace std;
using namespace ThePEG;

/** \ingroup Interfaces
 * 
 *  Some comment should be provided!
 */
class RootInterface {

 public:
    
    RootInterface(){}
    ~RootInterface(){}

    void init(string filename="Herwig++.root", string treename="tree");
    void finish();
    void fill(const string branchname, const double & value);
    void fillarray(const string branchname, const int &i, 
		   const double &px, const double &py, const double &pz, const double &E);
    void write();

 private:

    TFile *f;
    TTree *tree;
    TChain *chain;
    map <const string, TClonesArray*> array;
    map <const string, double> value;
};
 
}

#endif // ROOT_INTERFACE_H_
