#ifndef  __MgInterpolatorNested_h
#define  __MgInterpolatorNested_h

#include  "compvector.h"
#include  "mginterpolatorinterface.h"
#include  "multigridmeshinterface.h"
#include  "meshtransferinterface.h"
#include  <map>

/*-----------------------------------------*/


class MgInterpolatorNested : public MgInterpolatorInterface
{
protected:

  std::map<int,fixarray<2,int> >  zweier;
  std::map<int,fixarray<4,int> >  vierer;
  std::map<int,fixarray<8,int> >  achter;

  nvector<int>                c2f;

public:

  MgInterpolatorNested() {}
  ~MgInterpolatorNested() {}

  void init(const MeshTransferInterface* MT);
  
  void restrict_zero   (Vector&, const Vector&) const;
  void prolongate_add  (Vector&, const Vector&) const;
  void SolutionTransfer(Vector&, const Vector&) const;
  void Pi    (Vector& u) const;

};


#endif
