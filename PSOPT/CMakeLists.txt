cmake_minimum_required (VERSION 3.10)
project(virus_optimize LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/psopt-5.0.0/cmake)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(adolc REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(PSOPT REQUIRED)

pkg_check_modules(ipopt REQUIRED IMPORTED_TARGET ipopt)

add_executable(${PROJECT_NAME} $<$<BOOL:$<TARGET_EXISTS:PSOPT_SNOPT_interface>>:$<TARGET_OBJECTS:PSOPT_SNOPT_interface>> virus_optimize.cxx)
add_dependencies(${PROJECT_NAME} PSOPT)
target_link_libraries(${PROJECT_NAME} PUBLIC PSOPT PRIVATE PkgConfig::ipopt adolc Eigen3::Eigen)

if(${WITH_SNOPT_INTERFACE})
	enable_language(Fortran)
	target_link_libraries(${PROJECT_NAME} PRIVATE gfortran ${snopt7_LIBRARIES} PSOPT_SNOPT_interface)
endif()
