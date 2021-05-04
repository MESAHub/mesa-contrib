integer, parameter :: NUM_CZB_MESH_PARAM_SETS = 16

real(dp)       :: czb_mesh_f(NUM_CZB_MESH_PARAM_SETS)
real(dp)       :: czb_mesh_l(NUM_CZB_MESH_PARAM_SETS)
character(256) :: czb_mesh_zone_type(NUM_CZB_MESH_PARAM_SETS)
character(256) :: czb_mesh_zone_loc(NUM_CZB_MESH_PARAM_SETS)
character(256) :: czb_mesh_bdy_loc(NUM_CZB_MESH_PARAM_SETS)

namelist /czb_mesh_function/ &
     czb_mesh_zone_type, &
     czb_mesh_zone_loc, &
     czb_mesh_bdy_loc, &
     czb_mesh_f, &
     czb_mesh_l
  
