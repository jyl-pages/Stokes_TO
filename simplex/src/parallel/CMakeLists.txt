if(NOT USE_CUDA)
  return()
endif()
init_target(parallel)
use_ext_libs(LINK "cuda;")
use_src_libs(COMPILE "solver;")
setup_src_target()
