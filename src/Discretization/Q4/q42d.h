#ifndef __Q42d_h
#define __Q42d_h

#include "q4.h"

namespace Gascoigne
{

/**********************************************************/

  class Q42d : public Q4
  {
    protected:
      int GetPatchNumber(const Vertex2d& p0, Vertex2d& p) const;
      void VertexTransformation(const Vertex2d& p0, Vertex2d& p, int iq) const;

    public:
      Q42d();
      ~Q42d();

      std::string GetName() const {return "Q42d";}

      void BasicInit(const ParamFile* paramfile);
  };

/**********************************************************/

}

#endif
