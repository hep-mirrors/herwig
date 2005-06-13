// test program to check looptool linking against c++ code.

#include <iostream>
#include "Herwig++/Looptools/clooptools.h"

using namespace Herwig::Looptools;

int main() {
  ffini();
  std::cout << B0(1000.,50.,80.) << std::endl;
  std::cout << C0i(cc002,1000.,1000.,4500.,3.,5.,8.) << std::endl;
  ffexi();
  return 0;
}
