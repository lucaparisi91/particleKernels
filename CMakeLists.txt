cmake_minimum_required(VERSION 3.10)
# set the project name
project(particleKernels VERSION 3.0)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(COMPILER_FLAGS "")

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set (COMPILER_FLAGS "-mavx2 -mfma" )

else()
message( FATAL_ERROR "Compiler not supported (yet)." )
endif()


if ( ${CMAKE_BUILD_TYPE} MATCHES Debug)
  set(CMAKE_CXX_FLAGS_ADDITIONAL "-O0 -g -pg -Wfatal-errors ${COMPILER_FLAGS}")
  elseif (${CMAKE_BUILD_TYPE} MATCHES Release) 
  set(CMAKE_CXX_FLAGS_ADDITIONAL "-O3 -pg  -Wfatal-errors ${COMPILER_FLAGS}")
  else()
  message(FATAL_ERROR "Unrecognized build type: " ${CMAKE_BUILD_TYPE}  )

endif()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(MPI REQUIRED)

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS} ${CMAKE_CXX_FLAGS_ADDITIONAL} )
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

include_directories(SYSTEM ${MPI_INCLUDE_PATH})


add_subdirectory(${PROJECT_BINARY_DIR}/../external/googletest)


set(SOURCES  src/randomKernel.cpp src/cellNeighbourList.cpp src/timers.cpp
 src/lattice.cpp src/particleData.cpp src/utils.cpp)

 

 

add_library(main_lib STATIC ${SOURCES})
set_target_properties(main_lib PROPERTIES RULE_MESSAGES OFF)
set_target_properties( main_lib PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_ADDITIONAL} )
set_target_properties( main_lib PROPERTIES LINK_FLAGS "-pg" )


target_include_directories(main_lib PUBLIC
                           "${PROJECT_BINARY_DIR}/../src"
			   )

     	 
target_link_libraries (main_lib ${MPI_LIBRARIES} )


add_executable(main src/main.cpp)
target_link_libraries(main main_lib)
target_link_libraries(main )
set_target_properties( main PROPERTIES LINK_FLAGS "-pg" )
set_target_properties( main PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_ADDITIONAL} )


add_executable(tst tst/main.cpp tst/verlet-test.cpp tst/linkedCell.cpp)


set_target_properties( tst PROPERTIES LINK_FLAGS "-pg" )
set_target_properties( tst PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_ADDITIONAL} )
target_link_libraries(tst gtest main_lib ${MPI_LIBRARIES})		


#set(TEST_SOURCES  )
#add_executable(tst tst/main.cpp ${TEST_SOURCES})
#target_include_directories(tst PUBLIC
#                           "${PROJECT_BINARY_DIR}/../src"
#                           )
#target_link_libraries(tst PUBLIC ${PROJECT_BINARY_DIR}/../external gtest main_lib ${MPI_LIBRARIES})			 
#set_target_properties( tst PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_ADDITIONAL} )
#set_target_properties( tst PROPERTIES LINK_FLAGS "-pg" )

message("Compile Flags: " ${CMAKE_CXX_COMPILE_FLAGS} )