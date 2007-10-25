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
 *  This is a basic interface to write ROOT event files.
 * 
 */
class RootInterface {

 public:
    /** Standard constructor */
    RootInterface(){}
    /** Standard destructor */
    ~RootInterface(){}

    /** 
     * Initialization
     * @param filename is the desired event file name
     * @param treename is the name of the ROOT::TTree
    */
    void init(string filename="Herwig++.root", string treename="tree");

    /**
     * Method which has to be called at the end of the run to finalize the output.
     */
    void finish();

    /**
     * Fill a branch.
     * @param branchname is the name of the branch
     * @param value is the value that will be stored
     */
    void fill(const string branchname, const double & value);

    /**
     * Fill a TClonesArray of TLorentzVectors.
     * @param branchname is the name of the branch/TClonesArray
     * @param i is the counter
     * @param px is the x component of the TLorentzVector
     * @param py is the y component of the TLorentzVector
     * @param pz is the z component of the TLorentzVector
     * @param E is the E component of the TLorentzVector
     */
    void fillarray(const string branchname, const int &i, 
		   const double &px, const double &py, const double &pz, const double &E);

    /**
     * Method to write out one event. Has to be called after all branches 
     * or TClonesArrays have been filled.
     */
    void write();

 private:

    /** Pointer to the file */
    TFile *f;
    /** Pointer to the tree */
    TTree *tree;
    /** Pointer to a chain in case of multiple files */
    TChain *chain;
    /** map to store the TClonesArray's */
    map <const string, TClonesArray*> array;
    /** map to store the values for the regular branches */
    map <const string, double> value;
};
 
}

#endif // ROOT_INTERFACE_H_
