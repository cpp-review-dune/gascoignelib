#ifndef __integrationformulainterface_h
#define __integrationformulainterface_h

#include  "nvector.h"
#include  "vertex.h"

/*------------------------------------------------------------*/

namespace Gascoigne
{
  class IntegrationFormulaInterface
  {
    private:

    protected:

    public:
      IntegrationFormulaInterface() {}
      virtual ~IntegrationFormulaInterface() {}

      virtual int    n()      const=0;
      virtual double w(int k) const=0;

      virtual void xi(Vertex1d& v, int k) const {
        std::cerr << "\"IntegrationFormulaInterface::xi\" not written!" << std::endl;
        abort();
      }
      virtual void xi(Vertex2d& v, int k) const {
        std::cerr << "\"IntegrationFormulaInterface::xi\" not written!" << std::endl;
        abort();
      }
      virtual void xi(Vertex3d& v, int k) const {
        std::cerr << "\"IntegrationFormulaInterface::xi\" not written!" << std::endl;
        abort();
      }
  };
}

/*------------------------------------------------------------*/

#endif
