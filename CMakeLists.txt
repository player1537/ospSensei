## Copyright 2009-2020 Intel Corporation
## SPDX-License-Identifier: Apache-2.0

enable_language(C)


find_package(
  SENSEI
  REQUIRED
)


add_library(
  OSPRayStudioVisualization
SHARED
  src/OSPRayStudioVisualization.cpp
)


target_include_directories(
  OSPRayStudioVisualization
PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)


target_link_libraries(
  OSPRayStudioVisualization
PUBLIC
  sensei
  ospray_sg
)


install(
  TARGETS
    OSPRayStudioVisualization
  COMPONENT lib
)


add_executable(
  nbody
  app/nbody/nbody.cpp
)

target_link_libraries(
  nbody
PUBLIC
  OSPRayStudioVisualization
)


install(
  TARGETS
    nbody
  COMPONENT bin
)


option(BUILD_PLUGIN_SENSEI "Sensei plugin" OFF)

if (BUILD_PLUGIN_SENSEI)
  set( pluginName "ospray_studio_plugin_sensei")

  add_library(${pluginName} SHARED
    src/plugin_sensei.cpp
  )

  target_link_libraries(${pluginName} OSPRayStudioVisualization)

  target_link_libraries(${pluginName} ospray_sg)

  # Only link against imgui if needed (ie, pure file importers don't)
  target_link_libraries(${pluginName} imgui)

  target_include_directories(${pluginName}
    PRIVATE ${CMAKE_SOURCE_DIR}
  )

  install(TARGETS ${pluginName}
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
      COMPONENT lib
    # on Windows put the dlls into bin
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      COMPONENT lib
  )

endif()
