#ifndef  __VisuDataRhs_h
#define  __VisuDataRhs_h

#include  "visudata.h"
#include  "righthandsidedata.h"
#include  "meshinterface.h"

/////////////////////////////////////////////
///
///@brief
///  ... comments VisuDataRhs

///
///
/////////////////////////////////////////////

class VisuDataRhs : public VisuData
{
public:

private:

protected:

  const MeshInterface*               M;
  const RightHandSideData*  RHS;

public:

//
///  Constructor 
//

  VisuDataRhs(const RightHandSideData*  rhs, const MeshInterface* m);

  int    visucomp()     const;
  int    visun()        const;
  double visudata(int i,int c, const Vertex2d& v) const;
  double visudata(int i,int c, const Vertex3d& v) const;

};


#endif
