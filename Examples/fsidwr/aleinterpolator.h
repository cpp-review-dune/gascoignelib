/*----------------------------   aleinterpolator.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __aleinterpolator_H
#define __aleinterpolator_H
/*----------------------------   aleinterpolator.h     ---------------------------*/



#include  "mginterpolatornested.h"
#include  "alebasediscretization.h"

/*-----------------------------------------*/


namespace Gascoigne
{
  class AleInterpolator : public virtual MgInterpolatorNested
  {
  private:


    const HASHSET<int>      *__inodes;
    const HASHMAP<int,int>  *__fnodes, *__snodes;

  public:

  AleInterpolator() : MgInterpolatorNested() {}

    void SetInterface(const HASHSET<int>*     inodes,
		      const HASHMAP<int,int>* fnodes,
		      const HASHMAP<int,int>* snodes)
    {
      __inodes = inodes;
      __fnodes = fnodes;
      __snodes = snodes;
    }
    	      
    
    //    void BasicInit(const MeshTransferInterface* MT);
  
    void restrict_zero   (GlobalVector&, const GlobalVector&) const;
    void prolongate_add  (GlobalVector&, const GlobalVector&) const;
    /* void SolutionTransfer(GlobalVector&, const GlobalVector&) const; */
    /* void Pi    (GlobalVector& u) const; */

  };
}



/*----------------------------   aleinterpolator.h     ---------------------------*/
/* end of #ifndef __aleinterpolator_H */
#endif
/*----------------------------   aleinterpolator.h     ---------------------------*/
