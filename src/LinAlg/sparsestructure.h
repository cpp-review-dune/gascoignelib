#ifndef __sparsestructure_h
#define __sparsestructure_h

#include  <vector>
#include  "gascoigne.h"
#include  "sparsestructureinterface.h"

/*------------------------------------------------------------------------*/

namespace Gascoigne
{
class SparseStructure : public SparseStructureInterface
{
  protected:

    typedef  std::vector<IntSet>  Indices;    

    int       sntot;
    Indices   sindices;

  public:

  typedef  IntSet::iterator         iterator; 
  typedef  IntSet::const_iterator   const_iterator; 
  
  SparseStructure() : sntot(0),sindices(0) {}

  void clear()
  {
    sntot=0;
    sindices.clear();
  }
  

    int                n()              const { return sindices.size(); }
    int                ntotal()         const { return sntot; }

    const Indices&     indices()        const { return sindices; }
    const IntSet&      row(int i)       const 
      {
	assert((i>=0)&&(i<sindices.size()));
	return sindices[i];
      } 
    IntSet&             row(int i)
      {
	assert((i>=0)&&(i<sindices.size()));
	return  sindices[i];
      }
    IntSet::iterator            rowbegin(int i)
      {
	assert((i>=0)&&(i<sindices.size()));
	return  row(i).begin();
      }
    IntSet::iterator            rowend(int i)
      {
	assert((i>=0)&&(i<sindices.size()));
	return  row(i).end();
      }
    IntSet::const_iterator      rowbegin(int i)  const
      {
	assert((i>=0)&&(i<sindices.size()));
	return  row(i).begin();
      }
    IntSet::const_iterator      rowend(int i)    const
      {
	assert((i>=0)&&(i<sindices.size()));
	return  row(i).end();
      }
    int                rowsize(int i)   const
      {
	assert((i>=0)&&(i<sindices.size()));
	return  row(i).size();
      }

    friend std::ostream& operator<<(std::ostream &s, const SparseStructure& A);
    void statistics(std::ostream&) const;

    SparseStructure& operator=(const SparseStructure& B);
    void build_begin(int n);
    void build_clear(int i);
    void build_add(int i, int j)
      {
	row(i).insert(j);
      }
    template<class IT>
    void build_add(int i, IT lsf, IT lsl)
      {
	for(IT p=lsf;p!=lsl;p++) row(i).insert(*p);
      }
    template<class IT>
    void build_add(IT lsf, IT lsl)
      {
	for(IT p=lsf;p!=lsl;p++) build_add(*p,lsf,lsl);
      }
    template<class IT>
    void build_add(IT rf, IT rl, IT cf, IT cl)
      {
	for(IT p=rf;p!=rl;p++) build_add(*p,cf,cl);
      }
    void build_end();

    void hanging_node(int,int,int);
  
    void  enlarge(const SparseStructure&);
    void  enlarge_lu();
    void  enlarge_for_lu(const IntVector& perm);
};
}

#endif
