#ifndef __deletecell_h
#define __deletecell_h

#include  "gascoigne.h"

/*---------------------------------------------------*/

namespace Gascoigne
{
template <class C>
void delete_cells(const IntSet&, std::vector<C>&, const IntVector&, const IntVector&);
}

/*---------------------------------------------------*/

#endif
