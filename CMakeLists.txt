cmake_minimum_required(VERSION 2.8.10)

project(radprop)
set(libname "radprop")

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED)
find_package(GeographicLib REQUIRED)

#-- locate GDAL 
find_package(GDAL) 
#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
message("ROOT_VERSION is set to ${ROOT_VERSION}")
if( ${ROOT_VERSION} VERSION_GREATER "5.99/99")
    message("Using ROOT_VERSION 6")
    include(${ROOT_USE_FILE})
else()
     add_definitions(${ROOT_DEFINITIONS})	
endif()



include_directories(${PROJECT_SOURCE_DIR} ${GDAL_INCLUDE_DIRS} ${GEOGRAPHICLIB_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/src )

set (HEADER_FILES radprop/dem.h ) 
set (SOURCE_FILES src/dem.cc )  

set(DICTNAME G__${libname})


ROOT_GENERATE_DICTIONARY(${DICTNAME} ${HEADER_FILES} LINKDEF LinkDef.h)

#---Create a shared library with geneated dictionary
add_library(${libname} SHARED ${SOURCE_FILES} ${DICTNAME}.cxx)


if(GDAL_FOUND) 
  add_definitions( -DHAVE_GDAL )
  target_link_libraries(${libname} ${GDAL_LIBRARIES} ${ZLIB_LIBRARIES} ${ROOT_LIBRARIES} ${GEOGRAPHICLIB_LIBRARIES} )
else()
  target_link_libraries(${libname} ${ZLIB_LIBRARIES} ${ROOT_LIBRARIES} ${GEOGRAPHICLIB_LIBRARIES} )
endif() 




if( ${ROOT_VERSION} VERSION_GREATER "5.99.99")
    add_custom_target(${DICTNAME}.pcm DEPENDS ${DICTNAME})
endif()



set(CMAKE_CXX_FLAGS_DEFAULT 
  "-Os -g -Wall" 
  CACHE STRING "c++ Flags used during default builds" 
  FORCE ) 

set(CMAKE_C_FLAGS_DEFAULT 
  "-Os -g -Wall"
  CACHE STRING "c Flags used during default builds" 
  FORCE ) 

set(CMAKE_EXE_LINKER_FLAGS_DEFAULT 
  "-g"
  CACHE STRING "ld Flags used during default builds" 
  FORCE ) 

set(CMAKE_SHARED_LINKER_FLAGS_DEFAULT 
  "-g"
  CACHE STRING "ld Flags used during default builds" 
  FORCE ) 


mark_as_advanced ( CMAKE_CXX_FLAGS_DEFAULT CMAKE_C_FLAGS_DEFAULT CMAKE_EXE_LINKER_FLAGS_DEFAULT CMAKE_SHARED_LINKER_FLAGS_DEFAULT) 

if (NOT CMAKE_BUILD_TYPE) 
  set (CMAKE_BUILD_TYPE Default
    CACHE STRING "Choose tye type of build: None Debug Release RelWithDebInfo MinSizeRel Default"
    FORCE ) 
endif()







