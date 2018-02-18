project(CIDANE Fortran C CXX)
cmake_minimum_required (VERSION 2.8.2)

if(NOT DEFINED SEQAN_INCLUDE_PATH)
message( FATAL_ERROR "ERROR: \${SEQAN_INCLUDE_PATH} variable is not defined. Please define it.")    
#elseif(NOT DEFINED CPLEX_ROOT_DIR)
#message( FATAL_ERROR "ERROR: \${CPLEX_ROOT_DIR} variable is not defined. Please define it.")    
endif(NOT DEFINED SEQAN_INCLUDE_PATH)

if(NOT DEFINED CMAKE_SEQAN_MODULE_PATH)
get_filename_component(seqanParent ${SEQAN_INCLUDE_PATH} "DIRECTORY")
set(CMAKE_SEQAN_MODULE_PATH "${seqanParent}/util/cmake/")
message(STATUS "\${CMAKE_SEQAN_MODULE_PATH} variable not defined. Assuming FindSeqan.cmake to be located at: ${CMAKE_SEQAN_MODULE_PATH}")
endif(NOT DEFINED CMAKE_SEQAN_MODULE_PATH)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SEQAN_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}/cmake/Modules/")
set(local_dir ${CMAKE_HOME_DIRECTORY})

# set this here otherwise seqan find module will complain    
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEQAN_DISABLE_VERSION_CHECK -std=c++11")

set (SEQAN_FIND_DEPENDENCIES ZLIB BZip2)
find_package (SeqAn REQUIRED)

if (NOT ${SEQAN_HAS_ZLIB})
message(FATAL_ERROR "MISSING ZLIB")
endif(NOT ${SEQAN_HAS_ZLIB})

if (NOT ${SEQAN_HAS_BZip2})
message(FATAL_ERROR "MISSING ZLIB")
endif(NOT ${SEQAN_HAS_BZip2})

# Add flags for SeqAn (and its dependencies).
add_definitions (${SEQAN_DEFINITIONS})

#find_package(CPLEX REQUIRED)
find_package(GSL REQUIRED)

include_directories(
"${local_dir}/include/"
"${local_dir}/contrib/glmnet/"
${GSL_INCLUDE_DIR}
#${CPLEX_INCLUDE_DIRS}
${SEQAN_INCLUDE_DIRS}
)

include(source/sources.cmake)

set(LIBRARY_OUTPUT_PATH "./lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "./bin")

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -fexceptions -DIL_STD -Wall -Wextra -Wno-variadic-macros -pedantic")
set(CMAKE_Fortran_FLAGS "-c -ffixed-form")

#add_executable(CIDANE ./source/CIDANE.C ${sources})
#add_executable(exonRefine ./contrib/exonRefine/exonRefine.cpp)
add_executable(sigcount ./source/sigcount.cpp)

#target_link_libraries(CIDANE glmnet ${CPLEX_LIBRARIES} ${GSL_LIBRARIES} ${SEQAN_LIBRARIES} "-lm -pthread")
target_link_libraries(sigcount ${SEQAN_LIBRARIES})
#target_link_libraries(exonRefine ${SEQAN_LIBRARIES})

#### end
