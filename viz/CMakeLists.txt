set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(SpatialIndex 1.8.5 REQUIRED)
find_package(PNG 1.2.50 REQUIRED)
find_package(Boost 1.54 REQUIRED COMPONENTS
	system filesystem log program_options
)
find_package(Eigen3 3.2.2 REQUIRED)

add_executable(viz
	viz.cpp
	PNGWriter.hpp
	PNGWriter.cpp
	SMVStreamer.hpp
	SMVStreamer.cpp
	Reconstructer.hpp
	Reconstructer.cpp
)
target_include_directories(viz PUBLIC
  ${SPATIALINDEX_INCLUDE_DIR}
  ${PNG_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
)
target_link_libraries (viz LINK_PUBLIC
  ${SPATIALINDEX_LIBRARY}
  ${PNG_LIBRARIES}
  ${Boost_LIBRARIES}
)
add_definitions(
  ${PNG_DEFINITIONS}
)