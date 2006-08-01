#ifndef __Q43d_h
#define __Q43d_h

#include "q4.h"

namespace Gascoigne
{

/**********************************************************/

  class Q43d : public Q4
  {
    protected:
      int GetPatchNumber(const Vertex3d& p0, Vertex3d& p) const;
      void VertexTransformation(const Vertex3d& p0, Vertex3d& p, int iq) const;

    public:
      Q43d();
      ~Q43d();

      std::string GetName() const {return "Q43d";}

      void BasicInit(const ParamFile* paramfile);
  };

/**********************************************************/

}

#endif
