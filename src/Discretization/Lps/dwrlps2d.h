#ifndef __dwrlps2d_h
#define __dwrlps2d_h

#include "q1lps2d.h"

namespace Gascoigne
{
/*-------------------------------------------------*/

class DwrLps2d : public Q1Lps2d
{
 protected:

  const PatchMesh* GetPatchMesh() const {
    const PatchMesh* MP = dynamic_cast<const PatchMesh*>(GetMesh());
    assert(MP);
    return MP;
  }

 public:

  DwrLps2d() : Q1Lps2d() {}
  ~DwrLps2d() {}

  std::string GetName() const {return "DwrLps2d";}
  void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
};
}
/*-------------------------------------------------*/

#endif
