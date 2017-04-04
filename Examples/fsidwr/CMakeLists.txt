  
ADD_EXECUTABLE("newdwr" main.cc ale.cc aleloop.cc  chi.cc
pi-fsi.cc
aleinterpolator.cc
aleboundary.cc alesolver.cc alediscretization.cc dwr.cc aleadjoint.cc
ale-dual.cc alebase.cc 
aleblocks.cc
ale_slow.cc alemultilevelsolver.cc
 dwr_boundary.cc dwr_integrator.cc dwr_discretization.cc
)


ADD_EXECUTABLE("newdwr-3d" main_3d.cc ale.cc aleloop.cc  chi.cc
pi-fsi.cc
aleinterpolator.cc
aleboundary.cc alesolver.cc alediscretization.cc dwr.cc aleadjoint.cc
ale-dual.cc alebase.cc 
aleblocks.cc
ale_slow.cc alemultilevelsolver.cc
 dwr_boundary.cc dwr_integrator.cc dwr_discretization.cc
)



