#ifndef __visudata_h
#define __visudata_h

#include  "vertex.h"

/***************************************************************/

class VisuData
{
 protected:

 public:


  virtual ~VisuData(){}
  virtual int    visucomp()     const {return 0;}
  virtual int    visun()        const {return 0;}
  virtual double visudata(int i,int c) const {assert(0);}
  virtual double visudata(int i,int c, const Vertex2d& v) const {
    return visudata(i,c);
  }
  virtual double visudata(int i,int c, const Vertex3d& v) const {
    return visudata(i,c);
  }
};

/***************************************************************/

#endif
