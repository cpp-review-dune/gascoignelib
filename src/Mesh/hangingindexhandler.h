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

  std::map<int,fixarray<6,int> >   hnq4;
  std::map<int,fixarray<26,int> >  hnq4face;

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

  const std::map<int,fixarray<6,int> >*  GetQ4Structure()     const { return &hnq4;}
  const std::map<int,fixarray<26,int> >* GetQ4StructureFace() const { return &hnq4face;}
  std::map<int,fixarray<6,int> >*  GetQ4Structure()      { return &hnq4;}
  std::map<int,fixarray<26,int> >* GetQ4StructureFace()  { return &hnq4face;}
};
}

#endif
