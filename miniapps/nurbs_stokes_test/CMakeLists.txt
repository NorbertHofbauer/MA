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

add_mfem_miniapp(
  nurbs_stokes_test
  MAIN nurbs_stokes_test.cpp
  LIBRARIES mfem)

if (MFEM_ENABLE_TESTING)
  add_test(NAME nurbs_stokes_test_per_ser
    COMMAND $<TARGET_FILE:nurbs_stokes_test> -no-vis
    -m ${PROJECT_SOURCE_DIR}/mesh/pipe-nurbs-boundary-test_2.mesh)
endif()

