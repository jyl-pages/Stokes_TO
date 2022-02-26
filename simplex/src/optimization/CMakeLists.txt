init_target(optimization)
use_ext_libs(LINK "ipopt;nlopt;mma")
use_src_libs(COMPILE "reservoir;")
setup_src_target()
