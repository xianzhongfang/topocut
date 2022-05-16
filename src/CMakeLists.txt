aux_source_directory(. src)

find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(Eigen REQUIRED)
include_directories(${EIGEN_INCLUDE_DIR})


add_executable(topocut ${src})
target_link_libraries(topocut
  ${Boost_LIBRARIES}
  gmpxx
  gmp
  )
