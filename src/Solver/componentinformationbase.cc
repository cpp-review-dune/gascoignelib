#include  "componentinformationbase.h"
#include  "problemdescriptorinterface.h"
#include  "domainrighthandside.h"
#include  "compose_name.h" 


namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments ComponentInformation

  ///
  ///
  /////////////////////////////////////////////

std::string ComponentInformationBase::GetName() const {
  return "ComponentInformationBase";
}

const int ComponentInformationBase::GetNScalars     () const{
  ProblemDescriptorInterface* ppdi = GetProblemDescriptorInterface();
  assert( ppdi );

  int ncomps = -1;

  int ncomps_rhs=-1;
  // geht nicht :/
  // const Application* prhs   = ppdi->GetRightHandSide();
  // if( prhs!=NULL ){ ncomps = ncomps_rhs = prhs->GetNcomp(); }

  int ncomps_brhs=-1;
  const BoundaryRightHandSide* pbrhs   = ppdi->GetBoundaryRightHandSide();
  if( pbrhs!=NULL ){ ncomps = ncomps_brhs = pbrhs->GetNcomp(); }

  int ncomps_equation=-1;
  const Equation* peq   = ppdi->GetEquation();
  if( peq!=NULL ){ ncomps = ncomps_equation = peq->GetNcomp(); }

  int ncomps_boundaryequation=-1;
  const BoundaryEquation* pbeq   = ppdi->GetBoundaryEquation();
  if( pbeq!=NULL ){ ncomps = ncomps_boundaryequation = pbeq->GetNcomp(); }

  //int ncomps_dirichlet=-1;
  // doesn't have this method
  // const DirichletData* pdd   = ppdi->GetDirichletData();
  // if( pdd!=NULL ){ ncomps = ncomps_dirichlet = pdd->GetNcomp(); }

  int ncomps_initialcondition=-1;
  // doesn't have this method
  // const Application* pic   = ppdi->GetInitialCondition();
  // if( pic!=NULL ){ ncomps = ncomps_initialcondition = pic->GetNcomp(); }

  assert(0<ncomps);
  if(0<ncomps_rhs              ) assert(ncomps==ncomps_rhs             ) ;
  if(0<ncomps_brhs             ) assert(ncomps==ncomps_brhs            ) ;
  if(0<ncomps_equation         ) assert(ncomps==ncomps_equation        ) ;
  if(0<ncomps_boundaryequation ) assert(ncomps==ncomps_boundaryequation) ;
  if(0<ncomps_initialcondition ) assert(ncomps==ncomps_initialcondition) ;

  return ncomps;
}
void      ComponentInformationBase::GetScalarName   (int i, std::string& s_name) const{
  s_name="u";
  compose_name_without_dot(s_name,i); 
}
const int ComponentInformationBase::GetNVectors     () const{
  int ncomps = GetNcomp();
  if( ncomps<=2) return 0;
  return 1;
}
void      ComponentInformationBase::GetVectorName   (int i, std::string& s_name) const{
  s_name="v";
}
void      ComponentInformationBase::GetVectorIndices(int i, fixarray<3,int>& fa_vectorindices) const{
  if ( GetDimension() ==2) {
    fa_vectorindices[0] = 1;
    fa_vectorindices[1] = 2;
    fa_vectorindices[2] =-1;
  } else if ( GetDimension() ==3) {
    fa_vectorindices[0] = 1;
    fa_vectorindices[1] = 2;
    fa_vectorindices[2] = 3;
  } else {
    std::cerr << __FILE__ << " :bad dimension="<<GetDimension()<<"."<<std::endl;
    assert(0);
  }
}

}
