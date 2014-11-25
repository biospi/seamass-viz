set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(SpatialIndex 1.8.5 REQUIRED)
find_package(Boost 1.54 REQUIRED COMPONENTS
	system filesystem log program_options
)

add_executable(viz
	viz.cpp
	Reconstructer.hpp
	Reconstructer.cpp
	SMVStreamer.hpp
	SMVStreamer.cpp
)
target_include_directories(viz PUBLIC
  ${Boost_INCLUDE_DIRS}
  ${SPATIALINDEX_INCLUDE_DIR}
)
target_link_libraries (viz LINK_PUBLIC
  ${Boost_LIBRARIES}
  ${SPATIALINDEX_LIBRARY}
)
