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

  std::vector<Gascoigne::GlobalVector>  dat_node;

protected:


public:


//
////  Con(De)structor 
//

  LocalTimeLoop() : StdTimeLoop() {}
  void BasicInit(const Gascoigne::ParamFile* paramfile);

  void AddNodeVector(std::string filename);
  void NewMesh(const ProblemDescriptorInterface* PD);

  void init(std::string,int);
  void forward(std::string iname, int first, int last, const ProblemDescriptorInterface* PD);
  void backward(std::string iname, std::string name, int first, int last, const ProblemDescriptorInterface* PD);

};


#endif
