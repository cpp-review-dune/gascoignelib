/*----------------------------   pi-fsi.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __pi_fsi_H
#define __pi_fsi_H
/*----------------------------   pi-fsi.h     ---------------------------*/



#include  <map>
#include  "fixarray.h"
#include  "compvector.h"
#include  "gascoignemesh2d.h"
#include  "gascoignemesh3d.h"
#include  "alebasediscretization.h"

namespace Gascoigne
{
  /*-----------------------------------------*/

  class PiFSI
  {
  protected:
  
    std::map<int,fixarray<2,int> > edge;
    std::map<int,fixarray<4,int> > face;
    std::map<int,fixarray<8,int> > cell;

    void Init2d(const GascoigneMesh2d* MP);
    void Init3d(const GascoigneMesh3d* MP);
    const AleBaseDiscretization* __ALE;
    
  public:
  
  PiFSI(const AleBaseDiscretization* ALE):
    __ALE(ALE)
    {
    }
    

    void Init(const MeshInterface* MP);

    void vmult(CompVector<double>& y, const CompVector<double>& x, 
	       double s=1.) const;
    void vmultsolid(CompVector<double>& y, const CompVector<double>& x, 
		    double s=1.) const;
    void vmultfluid(CompVector<double>& y, const CompVector<double>& x, 
		    double s=1.) const;
    void vmultprimal(CompVector<double>& y, const CompVector<double>& x, 
		    double s=1.) const;
  };
}



/*----------------------------   pi-fsi.h     ---------------------------*/
/* end of #ifndef __pi_fsi_H */
#endif
/*----------------------------   pi-fsi.h     ---------------------------*/
