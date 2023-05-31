# Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.
#
# This file is part of the MFEM library. For more information and source code
# availability visit https://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

if (MFEM_USE_MPI)
  add_mfem_miniapp(
  nurbs_stokes_fsi_p
  MAIN nurbs_stokes_fsi.cpp
  nurbs_stokes_fsi_solver.cpp
  nurbs_stokes_fsi_solver.hpp
  viscosity_models.cpp
  viscosity_models.hpp
  interfacecoefficient.cpp
  interfacecoefficient.hpp
  LIBRARIES mfem)  
endif()

