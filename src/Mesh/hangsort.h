#ifndef __hangsort_h
#define __hangsort_h

/*---------------------------------------------------*/

class HangEdgeSort{
protected:
  const EdgeManager& HR;
public:
  HangEdgeSort(const EdgeManager& H) : HR(H) {}
  bool operator() (int i, int j) const
    {
      return !HR.EdgeIsHanging(i) && HR.EdgeIsHanging(j);
    }
};

/*---------------------------------------------------*/

class HangEdgeSort2{
protected:
  const EdgeManager& HR;
public:
  HangEdgeSort2(const EdgeManager& H) : HR(H) {}
  bool operator() (const Edge& i, const Edge& j) const
    {
      return !HR.EdgeIsHanging(i) && HR.EdgeIsHanging(j);
    }
};

#endif
