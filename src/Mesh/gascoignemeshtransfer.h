#ifndef  __GascoigneMeshTransfer_h
#define  __GascoigneMeshTransfer_h

#include  "meshtransferinterface.h"
#include  "fixarray.h"
#include  "nvector.h"
#include  <map>

/*-----------------------------------------*/

class GascoigneMeshTransfer : public MeshTransferInterface
{
protected:
  
  std::map<int,fixarray<2,int> >  zweier;
  std::map<int,fixarray<4,int> >  vierer;
  std::map<int,fixarray<8,int> >  achter;
  
  nvector<int>                c2f;
  std::map<int,int>             CellEiner;
  std::map<int,fixarray<4,int> >  CellVierer;
  std::map<int,fixarray<8,int> >  CellAchter;
  
public:
  
  const std::map<int,fixarray<2,int> >& GetZweier() const {return zweier;}
  const std::map<int,fixarray<4,int> >& GetVierer() const {return vierer;}
  const std::map<int,fixarray<8,int> >& GetAchter() const {return achter;}
  const nvector<int>&                   GetC2f()    const {return c2f;}
  
  std::map<int,fixarray<2,int> >& GetZweier() {return zweier;}
  std::map<int,fixarray<4,int> >& GetVierer() {return vierer;}
  std::map<int,fixarray<8,int> >& GetAchter() {return achter;}
  nvector<int>&                   GetC2f() {return c2f;}

  const std::map<int,int>             & GetCellEiner ()const  {return CellEiner;}
  const std::map<int,fixarray<4,int> >& GetCellVierer()const  {return CellVierer;}
  const std::map<int,fixarray<8,int> >& GetCellAchter()const  {return CellAchter;}

  std::map<int,int>             & GetCellEiner () {return CellEiner;}
  std::map<int,fixarray<4,int> >& GetCellVierer() {return CellVierer;}
  std::map<int,fixarray<8,int> >& GetCellAchter() {return CellAchter;}
  
  GascoigneMeshTransfer();
};

#endif
