#ifndef  __MgInterpolatorNested_h
#define  __MgInterpolatorNested_h

#include  "mginterpolatorinterface.h"
#include  "gascoigne.h"
#include  "meshtransferinterface.h"
#include  <map>

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
  
  void restrict_zero   (Gascoigne::GlobalVector&, const Gascoigne::GlobalVector&) const;
  void prolongate_add  (Gascoigne::GlobalVector&, const Gascoigne::GlobalVector&) const;
  void SolutionTransfer(Gascoigne::GlobalVector&, const Gascoigne::GlobalVector&) const;
  void Pi    (Gascoigne::GlobalVector& u) const;

};


#endif
