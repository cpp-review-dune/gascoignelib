#ifndef  __gnuplotdata_h
#define  __gnuplotdata_h

#include  <string>
#include "vertex.h"

/*-----------------------------------------*/

class GnuplotData
{
protected:

  std::string   plane;
  Vertex3d pos;

public:

  GnuplotData() {};
  GnuplotData(const std::string& s, const Vertex3d& pos);
  GnuplotData(const GnuplotData& GP);
  
  void   SetName(std::string& filename) const;
  bool   TestVertex(const Vertex2d& v) const;
  bool   TestVertex(const Vertex3d& v) const;
  double SetVertex (const Vertex2d& v) const;
  double SetVertex (const Vertex3d& v) const;
};

#endif
