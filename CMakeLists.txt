cmake_minimum_required(VERSION 2.8.10)

project(radprop)
set(libname "radprop")
set(CMAKE_CXX_STANDARD 11)

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



include_directories(${PROJECT_SOURCE_DIR} ${GDAL_INCLUDE_DIRS} ${GeographicLib_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/src )

set (HEADER_FILES radprop/dem.h radprop/coord.h radprop/prop.h) 
set (SOURCE_FILES src/dem.cc src/coord.cc src/prop.cc src/itwom3.0.cpp)  

set(DICTNAME G__${libname})


ROOT_GENERATE_DICTIONARY(${DICTNAME} ${HEADER_FILES} LINKDEF LinkDef.h)

#---Create a shared library with geneated dictionary
add_library(${libname} SHARED ${SOURCE_FILES} ${DICTNAME}.cxx)


if(GDAL_FOUND) 
  add_definitions( -DHAVE_GDAL )
  target_link_libraries(${libname} ${GDAL_LIBRARIES} ${ZLIB_LIBRARIES} ${ROOT_LIBRARIES} ${GeographicLib_LIBRARIES} )
else()
  target_link_libraries(${libname} ${ZLIB_LIBRARIES} ${ROOT_LIBRARIES} ${GeographicLib_LIBRARIES} )
endif() 

option ( ENABLE_OPENMP " Enable OpenMP for multicore usage") 
if(ENABLE_OPENMP) 
  FIND_PACKAGE( OpenMP REQUIRED)
  if(OPENMP_FOUND)
    message("OPENMP FOUND")
    add_definitions (-DENABLE_OPENMP) 
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  else() 
    message(FATAL_ERROR "you tried to use openmp but Cmake couldn't find compiler support") 
  endif() 
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







