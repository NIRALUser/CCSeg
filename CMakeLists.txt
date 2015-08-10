Project(CCSeg)
set(LOCAL_PROJECT_NAME CCSeg)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
set(BRAINS_CMAKE_HELPER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
include(${BRAINS_CMAKE_HELPER_DIR}/CMakeBuildMacros.cmake)
include(${BRAINS_CMAKE_HELPER_DIR}/PreventInSourceBuilds.cmake)
include(${BRAINS_CMAKE_HELPER_DIR}/CMakeCommonExternalDefinitions.cmake)

#-----------------------------------------------------------------------------
# Git protocole option
#-----------------------------------------------------------------------------
option(USE_GIT_PROTOCOL "If behind a firewall turn this off to use http
instead." ON)
set(git_protocol "git")
if(NOT USE_GIT_PROTOCOL)
  set(git_protocol "http")
endif()

#-----------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#External Projects
include(ExternalProject)
if(CMAKE_EXTRA_GENERATOR)
  set(gen "${CMAKE_EXTRA_GENERATOR} - ${CMAKE_GENERATOR}")
else()
  set(gen "${CMAKE_GENERATOR}")
endif()

#QT package
FIND_PACKAGE(Qt REQUIRED)
IF(QT_USE_FILE)
  INCLUDE_DIRECTORIES(${QT_INCLUDE_DIR})
  INCLUDE(${QT_USE_FILE})
ELSE(QT_USE_FILE)
   MESSAGE(FATAL_ERROR, "QT not found. Please set QT_DIR.")
ENDIF(QT_USE_FILE)

OPTION(COMPILE_CCSeg "Compile CCSeg." ON)
IF(COMPILE_CCSeg)
  # ITK must be searched for before SlicerExecutionModel as ITK_DEPEND needs to be initialized
  # in PACKAGE_NEEDS_ITK() before it is used in PACKAGE_NEEDS_SlicerExecutionModel()
  PACKAGE_NEEDS_ITK("${LOCAL_CMAKE_BUILD_OPTIONS}" "${gen}")
  PACKAGE_NEEDS_SlicerExecutionModel("${LOCAL_CMAKE_BUILD_OPTIONS}" "${gen}")
  PACKAGE_NEEDS_VTKWITHQT("${LOCAL_CMAKE_BUILD_OPTIONS}" "${gen}")
  set(proj CCSeg)
  ExternalProject_Add(${proj}
    DOWNLOAD_COMMAND ""
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${proj}
    BINARY_DIR ${proj}-build
    DEPENDS ${ITK_DEPEND} ${VTK_DEPEND} ${SlicerExecutionModel_DEPEND}
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      ${LOCAL_CMAKE_BUILD_OPTIONS}    
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      ${COMMON_EXTERNAL_PROJECT_ARGS}
      -DGenerateCLP_DIR:PATH=${GenerateCLP_DIR}
      -DITK_DIR:PATH=${ITK_DIR}
      -DVTK_DIR:PATH=${VTK_DIR}
    INSTALL_COMMAND ""
  )
  FORCE_BUILD_CHECK(${proj})
ENDIF(COMPILE_CCSeg)
