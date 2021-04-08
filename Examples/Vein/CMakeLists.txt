
ADD_EXECUTABLE("vein_3d" 
 Solver/alebase.cc    
 Solver/alediscretization.cc   
 Solver/solvers.cc 
 Solver/multilevelsolvers.cc 
 Solver/vankasolver.cc 
 Solver/vanka_matrix_vector.cc
 FSI/FSI_reduced/boundaryfsi.cc 
 FSI/FSI_reduced/fsi.cc   
 FSI/Deformation_Solid/def_solid.cc
 FSI/MeshMotion/meshmotion.cc 
 FSI/FSI_main/FSI_main.cc 
 FSI/FSI_main/FSI_main_dirichlet.cc 
 loop.cc
 main3d.cc  
 Mesh/finehierarchicalmesh3d.cc
 Mesh/localmeshagent.cc 
 Mesh/localhierarchicalmesh3d.cc
 Mesh/projectiononfinemeshagent.cc
 Fluid_Stat/fluid_stat.cc
 Fluid_Stat/boundary_fluid_stat.cc
 Fluid_Stat/fluid_stat_dirichlet.cc
 Solid_Euler/solideuler.cc
 Solid_Euler/boundarysolideuler.cc
 Solid_Euler/solid_euler_dirichlet.cc
)

SET_TARGET_PROPERTIES(
  vein_3d
  PROPERTIES
  COMPILE_DEFINITIONS "${GASCOIGNE_DEFINITIONS}"
  )

