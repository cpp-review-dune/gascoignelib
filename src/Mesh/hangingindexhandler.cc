#include  "hangingindexhandler.h"

using namespace std;

/*-----------------------------------------*/

HangingIndexHandler::HangingIndexHandler()
{
  hnq2.clear();
  hnq2face.clear();
}

/*-----------------------------------------*/

void HangingIndexHandler::CopyLevel2Nibble
(const HangingIndexHandler& Lhih,nvector<int>& Vg2l)
{
  hnq2.clear();
  hnq2face.clear();
      
  map<int,fixarray<3,int> >::const_iterator it3=Lhih.GetStructure()->begin();
  map<int,fixarray<9,int> >::const_iterator it9=Lhih.GetStructureFace()->begin();
  map<int,fixarray<3,int> >::const_iterator end3=Lhih.GetStructure()->end();
  map<int,fixarray<9,int> >::const_iterator end9=Lhih.GetStructureFace()->end();

  for (;it3!=end3;++it3)
    {
      int gf = it3->first;
      int lf = Vg2l[gf];
      if (lf<0) continue;
      fixarray<3,int> tmp;
      int gut=3;
      for (int i=0;i<3;++i)
	{
	  tmp[i]=Vg2l[it3->second[i]];
	  if (tmp[i]<0) --gut;
	  if (i<2) assert(tmp[i]>=0);
	}
      assert(gut==3);
      hnq2[lf]=tmp;
    }
  for (;it9!=end9;++it9)
    {
      int gf = it9->first;
      int lf = Vg2l[gf];
      if (lf<0) continue;
      fixarray<9,int> tmp;
      int gut=9;
      for (int i=0;i<9;++i)
	{
	  tmp[i]=Vg2l[it9->second[i]];
	  if (tmp[i]<0) --gut;
	}
      assert(gut==9);
      hnq2face[lf]=tmp;
    }
}
