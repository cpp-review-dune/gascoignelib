#ifndef __line_h
#define __line_h

#include  "cell.h"

namespace Gascoigne
{
class Line : public Cell<2>
{
 protected:

 public:

  Line(int l = 0, int f = -1) : Cell(l,f) {}

  int nnchild() { return 2;}
};
}

#endif
