#ifndef __cuthillmckee_h
#define __cuthillmckee_h

#include "columnstencil.h"


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

class CuthillMcKee 
{
    const ColumnStencil* S;
    
//     Vertex2d dir2d;
//     Vertex3d dir3d;
    int      dimension;
    std::vector<int> neighbors;

  public:
    CuthillMcKee (const StencilInterface *s);
    CuthillMcKee ();

    void Permutate    (nvector<int> &perm);
//     void Permutate    (nvector<int> &perm, const Vertex2d v);
//     void Permutate    (nvector<int> &perm, const Vertex3d v);
    
//     bool operator()(int i,int j) const;
};

#endif
