#ifndef __cuthillmckee_h
#define __cuthillmckee_h

#include "columnstencil.h"
#include "dynamicstencil.h"


/* CuthillMcKee
 *
 * macht nen CuthillMcKee fuer die ILU.
 * Bedeutet etwa, das immer Knoten mit
 * moeglichst wenig Nachbarn zuerst
 * genommen werden. Die Knoten einer Matrixzeile
 * stehen immer beieinander. Das soll die Bandbreite
 * klein halten.				
 *
 * So Gehts:
 *
 * CuthillMcKee cmc (UnstructuredStencilPointer);
 * cmc.Permutate (perm);   // reiner McKee
 * cmc.Permutate (perm,
 *                Vertex2d(1,0) ) // in x-Richtung
 * 
 */

namespace Gascoigne
{
class CuthillMcKee 
{
    const StencilInterface* S;
    const ColumnStencil*    CS;
    const DynamicStencil*   DS;
    
//     Vertex2d dir2d;
//     Vertex3d dir3d;
    int      dimension;
    std::vector<int> neighbors;

  public:
    CuthillMcKee (const StencilInterface *s);
    CuthillMcKee ();

    void Permutate    (IntVector &perm);
//     void Permutate    (IntVector &perm, const Vertex2d v);
//     void Permutate    (IntVector &perm, const Vertex3d v);
    
//     bool operator()(int i,int j) const;
    #ifdef __WITH_THREADS__
    void Permutate    (IntVector &perm, const IntVector &nodes_in_domain, const std::vector<std::vector<std::pair<int,int> > >& node2domain, int d);
    #endif
};
}

#endif
