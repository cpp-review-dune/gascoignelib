#ifndef  __MgInterpolatorNested_h
#define  __MgInterpolatorNested_h

#include  "mginterpolatorinterface.h"
#include  "gascoigne.h"
#include  "meshtransferinterface.h"
#include  <map>

using namespace Gascoigne;

/*-----------------------------------------*/


class MgInterpolatorNested : public virtual MgInterpolatorInterface
{
private:


  std::map<int,fixarray<2,int> >  zweier;
  std::map<int,fixarray<4,int> >  vierer;
  std::map<int,fixarray<8,int> >  achter;

  nvector<int>                c2f;

public:


  MgInterpolatorNested() : MgInterpolatorInterface() {}

  void init(const MeshTransferInterface* MT);
  
  void restrict_zero   (GlobalVector&, const GlobalVector&) const;
  void prolongate_add  (GlobalVector&, const GlobalVector&) const;
  void SolutionTransfer(GlobalVector&, const GlobalVector&) const;
  void Pi    (GlobalVector& u) const;

};


#endif
