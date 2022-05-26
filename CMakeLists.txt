## Copyright 2009-2020 Intel Corporation
## SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.22)
project(OSPSensei)

enable_language(C)


find_package(
  SENSEI
  REQUIRED
)


find_package(
  VTK
  REQUIRED
)


find_package(
  ospray
  REQUIRED
)


add_library(
  ospSensei
SHARED
  src/ospSensei/OSPRayVisualization.cpp
)


target_include_directories(
  ospSensei
PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)


target_link_libraries(
  ospSensei
PUBLIC
  ${VTK_LIBRARIES}
  sensei
  ospray::ospray
)

install(
  TARGETS
    ospSensei
  COMPONENT lib
)


add_executable(
  nbody
  app/nbody/nbody.cpp
)

target_link_libraries(
  nbody
  PUBLIC
    ospSensei
)

install(
  TARGETS
    nbody
  COMPONENT bin
)


add_executable(
  mandelbrot
  app/mandelbrot/mandelbrot.cpp
)

target_link_libraries(
  mandelbrot
  PUBLIC
    ospSensei
)

install(
  TARGETS
    mandelbrot
  COMPONENT bin
)
