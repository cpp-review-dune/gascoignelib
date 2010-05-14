#ifndef __NewmarkAlgorithm_h
#define __NewmarkAlgorithm_h

#include  "nonstationaryalgorithm.h"

/*----------------------------------------------------------------------------*/

class NewmarkAlgorithm : public Gascoigne::NonstationaryAlgorithm
{
public:
  
  void Run(const std::string& problemlabel);
};

#endif
