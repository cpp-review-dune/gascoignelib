#ifndef  __MgInterpolatorInterface_h
#define  __MgInterpolatorInterface_h

#include  "gascoigne.h"

/*--------------------------------------------------------*/

namespace Gascoigne
{
  class MgInterpolatorInterface
  {
    private:

    protected:

    public:
      MgInterpolatorInterface() {}
      virtual ~MgInterpolatorInterface() {}

      virtual void restrict_zero   (GlobalVector&, const GlobalVector&) const=0;
      virtual void prolongate_add  (GlobalVector&, const GlobalVector&) const=0;
      virtual void SolutionTransfer(GlobalVector&, const GlobalVector&) const=0;
      virtual void SolutionTransferUp(GlobalVector& ul, const GlobalVector& uL) const {
        prolongate_add(ul,uL); 
      }
      virtual void Pi(GlobalVector& u) const {
        std::cerr << "\"MgInterpolatorInterface::Pi\" not written!" << std::endl;
        abort();
      }
  };
}

#endif
