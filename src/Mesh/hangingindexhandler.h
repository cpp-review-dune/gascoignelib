#ifndef  __HangingIndexHandler_h
#define  __HangingIndexHandler_h

#include  "fixarray.h"
#include  <map>
#include  "gascoigne.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class HangingIndexHandler
{
 protected:
  
  typedef  fixarray<2,int>  IntVector2;
  typedef  fixarray<4,int>  IntVector4;
  
  std::map<int,fixarray<3,int> >  hnq2;
  std::map<int,fixarray<9,int> >  hnq2face;

 public:

  HangingIndexHandler();
  
  void Equal(const std::map<int,fixarray<3,int> >& h2,
	     const std::map<int,fixarray<9,int> >& h2f) 
    {
      hnq2=h2; 
      hnq2face=h2f;
    }

  void CopyLevel2Nibble
  (const HangingIndexHandler& Lhih,IntVector& Vg2l);

  // zugriff

  const std::map<int,fixarray<3,int> >* GetStructure()     const { return &hnq2;}
  const std::map<int,fixarray<9,int> >* GetStructureFace() const { return &hnq2face;}
  std::map<int,fixarray<3,int> >* GetStructure()      { return &hnq2;}
  std::map<int,fixarray<9,int> >* GetStructureFace()  { return &hnq2face;}
};
}

#endif
