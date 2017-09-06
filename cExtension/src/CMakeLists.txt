# ================================================================================
# 
# This file is part of the MagBottleC project.
# MagBottleC : libMagBottleC
# 
# File: CMakeLists.txt
# 
# Original Author(s) of this File:
# 	Tom Gade, 2017
# 
# Author(s) of Significant Updates/Modifications to the File:
# 	...
# 
# -----------------------------------------------------------------------------------

project(libMagBottleC)

set(mbc_bfield_cpp
  ${mbc_src_dir}/bfield/BField_c.cpp
  ${mbc_src_dir}/bfield/WireCoil_c.cpp
)

set(mbc_bfield_h
  ${mbc_inc_dir}/bfield/BField_c.h
  ${mbc_inc_dir}/bfield/BObject_c.h
  ${mbc_inc_dir}/bfield/WireCoil_c.h
)

set(mbc_integrate_cpp
  ${mbc_src_dir}/integrate/gauss_legendre.c
)

set(mbc_integrate_h
  ${mbc_inc_dir}/integrate/gauss_legendre.h
  ${mbc_inc_dir}/integrate/numerical_int.h
)

set(mbc_particle_cpp
  ${mbc_src_dir}/particle/Particle_c.cpp
)

set(mbc_particle_h
  ${mbc_inc_dir}/particle/Particle_c.h
)

set(mbc_tools_cpp
  ${mbc_src_dir}/tools/constsandtypes.cpp
  ${mbc_src_dir}/tools/numericalTechniques.cpp
  ${mbc_src_dir}/tools/vectortools.cpp
)

set(mbc_tools_h
  ${mbc_inc_dir}/tools/constsandtypes.h
  ${mbc_inc_dir}/tools/vectortools.h
)

set(mbc_vis_cpp
  ${mbc_src_dir}/vis/ItemDraw_c.cpp
  ${mbc_src_dir}/vis/MinVR.cpp
)

set(mbc_vis_h
  ${mbc_inc_dir}/vis/ItemDraw_c.h
  ${mbc_inc_dir}/vis/MinVR.h
)

set(mbc_api_cpp
	${mbc_src_dir}/api.cpp
)

set(mbc_api_h
	${mbc_inc_dir}/api.h
)

set(mbc_sources 
  ${mbc_bfield_cpp} ${mbc_bfield_h}
  ${mbc_integrate_cpp} ${mbc_integrate_h}
  ${mbc_particle_cpp} ${mbc_particle_h}
  ${mbc_tools_cpp} ${mbc_tools_h}
  ${mbc_api_cpp} ${mbc_api_h}
  #${mbc_vis_cpp} ${mbc_vis_h}
)

source_group("BField\\Header Files" FILES ${mbc_bfield_h})
source_group("BField\\Source Files" FILES ${mbc_bfield_cpp})
source_group("Integrate\\Header Files" FILES ${mbc_integrate_h})
source_group("Integrate\\Source Files" FILES ${mbc_integrate_cpp})
source_group("Particle\\Header Files" FILES ${mbc_particle_h})
source_group("Particle\\Source Files" FILES ${mbc_particle_cpp})
source_group("Tools\\Header Files" FILES ${mbc_tools_h})
source_group("Tools\\Source Files" FILES ${mbc_tools_cpp})
source_group("" FILES ${mbc_api_cpp} ${mbc_api_h})
#source_group("Vis\\Header Files" FILES ${mbc_vis_h})
#source_group("Vis\\Source Files" FILES ${mbc_vis_cpp})

add_library(MagBottleC SHARED ${mbc_sources})

install(TARGETS MagBottleC
        LIBRARY DESTINATION ${mbc_base_dir}/../lib
	ARCHIVE DESTINATION ${mbc_base_dir}/build/lib)

#install(FILES ${mbc_bfield_h} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/bfield)
#install(FILES ${mbc_integrate_h} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/integrate)
#install(FILES ${mbc_particle_h} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/particle)
#install(FILES ${mbc_tools_h} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/tools)
#install(FILES ${mbc_api_h} DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
#install(FILES ${mbc_vis_h} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/vis)

