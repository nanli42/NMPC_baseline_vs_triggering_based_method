cmake_minimum_required(VERSION 2.8.3)
project(iccr_demo)

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

## Build talker and listener
include_directories(
  /usr/local/include
  external/qpoases/
  external/qpoases/INCLUDE/
  external/qpoases/SRC/
  external/mpcc/spline/
  include/
)

## Declare a C++ library
# find ACADO generated code
file(GLOB ACADO_GENERATED_SRC_N15
  "src/acado_generated_code_N15/*.cpp"
  "src/acado_generated_code_N15/*.c"
  "external/qpoases/SRC/*.cpp"
  "external/mpcc/spline/*.cpp"
)
file(GLOB ACADO_GENERATED_SRC_N30
  "src/acado_generated_code_N30/*.cpp"
  "src/acado_generated_code_N30/*.c"
  "external/qpoases/SRC/*.cpp"
  "external/mpcc/spline/*.cpp"
)
add_library(acado_generated_rti_N15 STATIC ${ACADO_GENERATED_SRC_N15})
add_library(acado_generated_rti_N30 STATIC ${ACADO_GENERATED_SRC_N30})
target_include_directories(acado_generated_rti_N15 PUBLIC
  src/acado_generated_code_N15
)
target_include_directories(acado_generated_rti_N30 PUBLIC
  src/acado_generated_code_N30
)
file(GLOB SRC_track1_baseline
  "src/main-track1-baseline.cpp"
  "src/track1.cpp"
  "src/types.cpp"
)
file(GLOB SRC_track2_baseline
  "src/main-track2-baseline.cpp"
  "src/track2.cpp"
  "src/types.cpp"
)


link_directories(/usr/local/lib)

add_executable(main_track1_N15_baseline ${SRC_track1_baseline})
target_link_libraries(main_track1_N15_baseline acado_generated_rti_N15)
add_executable(main_track1_N30_baseline ${SRC_track1_baseline})
target_link_libraries(main_track1_N30_baseline acado_generated_rti_N30)
add_executable(main_track2_N15_baseline ${SRC_track2_baseline})
target_link_libraries(main_track2_N15_baseline acado_generated_rti_N15)
add_executable(main_track2_N30_baseline ${SRC_track2_baseline})
target_link_libraries(main_track2_N30_baseline acado_generated_rti_N30)

file(GLOB SRC_track1_triggering
  "src/main-track1-triggering.cpp"
  "src/track1.cpp"
  "src/types.cpp"
)
file(GLOB SRC_track2_triggering
  "src/main-track2-triggering.cpp"
  "src/track2.cpp"
  "src/types.cpp"
)

add_executable(main_track1_triggering ${SRC_track1_triggering})
target_link_libraries(main_track1_triggering acado_generated_rti_N30)
add_executable(main_track2_triggering ${SRC_track2_triggering})
target_link_libraries(main_track2_triggering acado_generated_rti_N30)
