#include  "localsolver.h"
#include  "localmeshinterpretor.h"

using namespace std;

/* ----------------------------------------- */

LocalSolver::LocalSolver() : StdSolver()
{
}

/* ----------------------------------------- */

MeshInterpretorInterface* LocalSolver::NewMeshInterpretor(int dimension, const std::string& discname)
{
  return new LocalMeshInterpretor;
}
