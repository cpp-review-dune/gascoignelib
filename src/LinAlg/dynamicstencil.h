#ifndef __dynamicstencil_h
#define __dynamicstencil_h

#include  <list>
#include  <vector>
#include  <cassert>

#include  "stencilinterface.h"


namespace Gascoigne
{

    class DynamicStencil : public StencilInterface
    {
    protected:
	typedef std::list<int>::const_iterator    const_citerator;
	typedef std::list<int>::iterator          citerator;
	
	// the column for every entry in the row

    public:
	std::vector<std::list<int> > cols; 
	
	int n() const { return cols.size(); }

	const_citerator cstart(int i) const { assert(i<n()); return cols[i].begin(); }
	const_citerator cstop (int i) const { assert(i<n()); return cols[i].end(); }
	citerator       cstart(int i)       { assert(i<n()); return cols[i].begin(); }
	citerator       cstop (int i)       { assert(i<n()); return cols[i].end(); }

    };
}


#endif
