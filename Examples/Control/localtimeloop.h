#ifndef  __LocalTimeLoop_h
#define  __LocalTimeLoop_h

#include  "stdtimeloop.h"


/////////////////////////////////////////////
////
////@brief
////  ... comments LocalTimeLoop

////
////
/////////////////////////////////////////////



class LocalTimeLoop : public StdTimeLoop
{
private:

  vector<GlobalVector>  dat_node;

protected:


public:


//
////  Con(De)structor 
//

  LocalTimeLoop() : StdTimeLoop() {}
  void BasicInit(const ParamFile* paramfile, const ProblemDescriptorInterface* PD);

  void AddNodeVector(string filename);
  void NewMesh();

  void init(string,int);
  void forward(string iname, int first, int last);
  void backward(string iname, string name, int first, int last);

};


#endif
