#ifndef __visudata_h
#define __visudata_h

#include  "vertex.h"

/***************************************************************/

namespace Gascoigne
{
  class VisuData
  {
    private:
      
    protected:

    public:
      VisuData() {}
      virtual ~VisuData(){}

      virtual int    visucomp()     const {
        return 0;
      }
      virtual int    visun()        const {
        return 0;
      }
      virtual double visudata(int i,int c) const {
        std::cerr << "\"VisuData::visudata\" not written!" << std::endl;
        abort();
      }
      virtual double visudata2(int i,int c, const Vertex2d& v) const {
        return visudata(i,c);
      }
      virtual double visudata2(int i,int c, const Vertex3d& v) const {
        return visudata(i,c);
      }
  };
}

/***************************************************************/

#endif
