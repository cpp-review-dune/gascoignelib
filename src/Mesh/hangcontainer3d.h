#ifndef __hangcontainer3d_h
#define __hangcontainer3d_h

#include "hangcontainer2d.h"

/*********************************************************************/

namespace Gascoigne
{
class HangContainer3d : public HangContainer2d
{
 public:
  
  typedef fixarray<4,int>    FaceVector;

  HangList<4>   FaceToBeDeleted, FaceToBeCreated;
  HangList<4>   FaceNotAnyMore;
  HangList<4>&  FaceHanging;

 public:

  HangContainer3d(HangList<2>& lh2, HangList<4>& lh3);

  const HangList<4>& FaceCreating() const { return FaceToBeCreated;}
  const HangList<4>& FaceDeleting() const { return FaceToBeDeleted;}
  const HangList<4>& FaceNotMore()  const { return FaceNotAnyMore;}
        HangList<4>& FaceNotMore()        { return FaceNotAnyMore;}

  void make_consistent();
  void output() const;

  int nFaceVertexesToBeDeleted() const { return FaceToBeDeleted.size(); }

  int nDel()  const { return VertexToBeDeleted.size()+FaceToBeDeleted.size(); }
  int nNew()  const { return VertexToBeCreated.size()+FaceToBeCreated.size(); }

  void load_elimination (IntVector&) const;
  int  vertex_index     (const EdgeVector& v) const
    {
      return HangContainer2d::vertex_index(v);
    }
  int  vertex_index(const FaceVector&) const;

  void update_olds(IntVector&, const IntVector&);
  void update_news(const IntVector&, int);

  void face_coarse(const FaceVector&, int, int);
  void face_refine(const FaceVector&, int);	
			      	   	
  void line_coarse(EdgeVector&, int, int);
  void line_refine(EdgeVector&, int, const HangList<2>& oldhangs);

  void NeighbourSwapper();
  void clear_hanging_lines();
  void build_hanging_lines(const HangList<2>& oldhangs);

  /* for boundary lines */
  bool ToBeDeleted(const FaceVector& v) const;
  bool ToBeCreated(const FaceVector& v) const;
};
}

/*********************************************************************/

#endif
