#include "Herwig++/Interfaces/AmegicInterface.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/EventRecord/Particle.h" 
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/Repository.h"
#include "Run_Parameter.H"
#include "Message.H"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;
using namespace ThePEG;
//using namespace AMEGIC;

// Using external variables from AMEGIC
AMEGIC::Run_Parameter rpa;
AMEGIC::Message AMEGIC::msg;

ParticleVector AmegicInterface::OneEvent(bool anti) 
{
   ParticleVector rval;
   vec4d *p;
   double weight;

   if(!initialized) initializeAmegic();
   weight = process->One_Event(p);
   if(HERWIG_DEBUG_LEVEL >= HwDebug::full) {
      for(int i = 0; i<inParticles.size() + outParticles.size(); i++) {
         generator()->log() << "The point " << p[i] << " has weight " << weight << endl;
      }
   }      

   for(int i = inParticles.size(); i<inParticles.size() + outParticles.size(); i++) {
      if(anti && outParticles[i-inParticles.size()]->CC()) {
         rval.push_back(outParticles[i-inParticles.size()]->CC()->produceParticle(
                                     LorentzMomentum(p[i][1]*GeV, p[i][2]*GeV, p[i][3]*GeV, p[i][0]*GeV)));
      } else {
         rval.push_back(outParticles[i-inParticles.size()]->produceParticle(
                                     LorentzMomentum(p[i][1]*GeV, p[i][2]*GeV, p[i][3]*GeV, p[i][0]*GeV)));
      }
   }
   return rval;
}

AmegicInterface::~AmegicInterface() 
{
//   if(process != NULL) delete process;
}

ClassDescription<AmegicInterface> AmegicInterface::initAmegicInterface;  

AmegicInterface & AmegicInterface::operator=(const AmegicInterface &a) { return *this; }

IBPtr AmegicInterface::clone() const { return new_ptr(*this); }
IBPtr AmegicInterface::fullclone() const { return clone(); }

void AmegicInterface::persistentOutput(PersistentOStream &os) const {
   os << setupDirectory << processString << inParticles << outParticles;
}

void AmegicInterface::persistentInput(PersistentIStream &is, int i) {
   is >> setupDirectory >> processString >> inParticles >> outParticles;
}
 
void AmegicInterface::Init() {
   static ClassDocumentation<AmegicInterface> documentation
     ("Class to interface the AMEGIC system with HERWIG++");
}

void AmegicInterface::readSetup(istream &is) throw(SetupException) {
   string temp;
   is >> temp;
   setProcess(temp);
 
   is >> temp;
   setupDirectory = temp;
}

void AmegicInterface::setProcess(string theProcessString) {
   int next;
   processString = theProcessString;
 
   next = theProcessString.find("->");
   if(next == -1) cerr << "Error in input decay string!!!\n";
   inParticles = createParticles(theProcessString.substr(0, next)+';');
   outParticles = createParticles(theProcessString.substr(next+2,theProcessString.find(';')));
}

PDVector AmegicInterface::createParticles(string list) {
   PDVector rval; 
   string temp;
   int next;
   do {
      next = min(list.find(','), list.find(';'));
      temp = list.substr(0,next);
      if(Repository::findParticle(temp)) 
         rval.push_back(Repository::findParticle(temp));
      else cerr << "Error in children input stream " << temp << endl;
      if(list.size() > next+1) list = list.substr(next+1);
      else break;
   } while(list[0] != ';' && list.size());
   return rval;
}

void AmegicInterface::doupdate() throw(UpdateException) {}
void AmegicInterface::doinit()  throw (InitException) {}
void AmegicInterface::dofinish() {}

void AmegicInterface::initializeAmegic() {
   Flavour *fin, *fout;
   uint i;
   AMEGIC::msg.SetFile();
   fin = new Flavour[inParticles.size()];
   fout = new Flavour[outParticles.size()];
   for(i = 0; i<inParticles.size(); i++) {
      fin[i] = AmegicInterface::convertParticle(inParticles[i]);
   }
   for(i = 0; i<outParticles.size(); i++) {
      fout[i] = AmegicInterface::convertParticle(outParticles[i]);
   }

   rpa.Init(setupDirectory);
   AMEGIC::particle_init(setupDirectory);
   process = new Amegic(setupDirectory, inParticles.size(), outParticles.size(), fin, fout);

   process->Prepare_Events();   
   initialized = true;
   process->Run();
}

Flavour AmegicInterface::convertParticle(PDPtr p) {
#define cmp(a,b,c) case ParticleID::a : return Flavour(kf::c); case ParticleID::b : return Flavour(kf::c).bar()
#define cmp2(a,c) case ParticleID::a : return Flavour(kf::c);
   switch (p->id()) {
      cmp(d,dbar,d);
      cmp(u,ubar,u);
      cmp(s,sbar,s);
      cmp(c,cbar,c);
      cmp(b,bbar,b);
      cmp(t,tbar,t);
      cmp(eminus,eplus,e);
      cmp(nu_e, nu_ebar, nue);
      cmp(muminus,muplus,mu);
      cmp(nu_mu, nu_mubar, numu);
      cmp(tauminus, tauplus, tau);
      cmp(nu_tau, nu_taubar, nutau);
      cmp2(g,gluon);
      cmp2(gamma,photon);
      cmp2(Z0,Z);
      cmp(Wminus,Wplus,W);
      cmp2(h0,h);
      cmp2(H0,H0);
      cmp2(A0,A0);
      cmp(SUSY_d_L,SUSY_d_Lbar,sDownL);
      cmp(SUSY_u_L,SUSY_u_Lbar,sUpL);
      cmp(SUSY_s_L,SUSY_s_Lbar,sStrangeL);
      cmp(SUSY_c_L,SUSY_c_Lbar,sCharmL); 
      cmp(SUSY_b_1,SUSY_b_1bar,sBottomL);
      cmp(SUSY_t_1,SUSY_t_1bar,sTopL);
      cmp(SUSY_e_Lminus, SUSY_e_Lplus, sElectronL);
      cmp(SUSY_nu_eL, SUSY_nu_eLbar,sNu1);
      cmp(SUSY_mu_Lminus, SUSY_mu_Lplus, sMuL);
      cmp(SUSY_nu_muL, SUSY_nu_muLbar, sNu2);
      cmp(SUSY_tau_1minus, SUSY_tau_1plus, sTauL);
      cmp(SUSY_nu_tauL, SUSY_nu_tauLbar, sNu3);
      cmp2(SUSY_g, Gluino);
      cmp(SUSY_chi_1plus, SUSY_chi_1minus, Chargino1);
      cmp(SUSY_chi_2plus, SUSY_chi_2minus, Chargino2);
      cmp2(SUSY_chi_10, Neutralino1);
      cmp2(SUSY_chi_20, Neutralino2);
      cmp2(SUSY_chi_30, Neutralino3);
      cmp2(SUSY_chi_40, Neutralino4);
      cmp(SUSY_d_R, SUSY_d_Rbar, sDownR);
      cmp(SUSY_u_R, SUSY_u_Rbar, sUpR);
      cmp(SUSY_s_R, SUSY_s_Rbar, sStrangeR);
      cmp(SUSY_c_R, SUSY_c_Rbar, sCharmR);
      cmp(SUSY_b_2, SUSY_b_2bar, sBottomR);
      cmp(SUSY_t_2, SUSY_t_2bar, sTopR);
      cmp(SUSY_e_Rminus, SUSY_e_Rplus, sElectronR);
      cmp(SUSY_mu_Rminus, SUSY_mu_Rplus, sMuR); 
      cmp(SUSY_tau_2minus, SUSY_tau_2plus, sTauR);
      cmp2(undefined,none);
   }
#undef cmp
#undef cmp2
}
