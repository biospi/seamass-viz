# Find Spatialindex
# ~~~~~~~~
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#
# Once run this will define:
#
# SPATIALINDEX_FOUND       = system has Spatialindex lib
# SPATIALINDEX_LIBRARY     = full path to the Spatialindex library
# SPATIALINDEX_INCLUDE_DIR = where to find headers
#


FIND_PATH(SPATIALINDEX_INCLUDE_DIR NAMES SpatialIndex.h PATHS
  /usr/include
  /usr/local/include
  "$ENV{LIB_DIR}/include"
  "$ENV{INCLUDE}"
  "$ENV{OSGEO4W_ROOT}/include"
  "${OSGEO4W_ROOT}/include"
  PATH_SUFFIXES spatialindex
  )

FIND_LIBRARY(SPATIALINDEX_LIBRARY NAMES spatialindex_i spatialindex spatialindex-32 spatialindex-64 PATHS
  /usr/lib
  /usr/local/lib
  "$ENV{LIB_DIR}/lib"
  "$ENV{LIB}/lib"
  "$ENV{OSGEO4W_ROOT}/lib"
  "${OSGEO4W_ROOT}/lib"
  )

IF (SPATIALINDEX_INCLUDE_DIR AND SPATIALINDEX_LIBRARY)
  SET(SPATIALINDEX_FOUND TRUE)
ENDIF (SPATIALINDEX_INCLUDE_DIR AND SPATIALINDEX_LIBRARY)

IF (SPATIALINDEX_FOUND)
  IF (NOT SPATIALINDEX_FIND_QUIETLY)
    MESSAGE(STATUS "Found libSpatialIndex: ${SPATIALINDEX_LIBRARY}")
  ENDIF (NOT SPATIALINDEX_FIND_QUIETLY)
ELSE (SPATIALINDEX_FOUND)
  IF (SPATIALINDEX_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Could not find libSpatialIndex")
  ENDIF (SPATIALINDEX_FIND_REQUIRED)
ENDIF (SPATIALINDEX_FOUND)