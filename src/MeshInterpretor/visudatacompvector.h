#ifndef __visudatacompvector_h
#define __visudatacompvector_h

#include  "visudata.h"
#include  "compvector.h"

/***************************************************************/

// alle ncomps gleich grosss !!!!!!!!!

class VisuDataCompVector : public VisuData
{
 protected:

  typedef const CompVector<double>* VectorPointer;
  std::vector<VectorPointer> vvp;
  std::vector<int> isizes;

  std::pair<int,int> GetIndex(int c) const;


 public:

  virtual ~VisuDataCompVector(){}
  VisuDataCompVector();
  VisuDataCompVector(const CompVector<double>& v);

  void AddVector(const CompVector<double>& v);
  virtual int    visucomp()     const;
  int    visun()        const;
  virtual double visudata(int i,int c) const;
};

/***************************************************************/

#endif
