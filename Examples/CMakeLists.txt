SUBDIRS(Laplace2D Laplace3D NavierStokes2D NavierStokes3D Heat2D ConvectionDiffusion2D Control)
LINK_LIBRARIES(Solver Problems MeshInterpretor LinAlg Mesh Common ${EXT_LIBS})

