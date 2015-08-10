
#-----------------------------------------------------------------------------
## The FORCE_BUILD_CHECK macro adds a forecebuild step that will cause the
## external project build process to be checked for updates each time
## a dependent project is built.  It MUST be called AFTER the ExternalProject_Add
## step for the project that you want to force building on.
macro(FORCE_BUILD_CHECK  proj)
    ExternalProject_Add_Step(${proj} forcebuild
      COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BUILD_DIR}/${proj}-prefix/src/${proj}-stamp/${proj}-build
      DEPENDEES configure
      DEPENDERS build
      ALWAYS 1
    )
endmacro()

#-----------------------------------------------------------------------------
## empty until ITK is brought into here as an ExternalProject
macro(PACKAGE_NEEDS_ITK LOCAL_CMAKE_BUILD_OPTIONS gen)
  set(packageToCheck ITK)
  OPTION(OPT_USE_SYSTEM_${packageToCheck} "Use the system's ${packageToCheck} library." OFF)
  #  MARK_AS_ADVANCED(OPT_USE_SYSTEM_${packageToCheck})
  if(OPT_USE_SYSTEM_ITK)
    find_package(ITK REQUIRED)
    include(${ITK_USE_FILE})
    set(ITK_DEPEND "") ## Set the external depandancy for ITK
  else()
    set(proj Insight)
    set(${proj}_REPOSITORY ${git_protocol}://itk.org/ITK.git)
    set(${proj}_GIT_TAG 042875b0246b0b2da1c1c58bf2f606b9590f6979 )
    set(ITK_VERSION_ID ITK-4.8)
    ExternalProject_Add(${proj}
      GIT_REPOSITORY ${${proj}_REPOSITORY}
      GIT_TAG ${${proj}_GIT_TAG}
      UPDATE_COMMAND ""
      SOURCE_DIR ${proj}
      BINARY_DIR ${proj}-build
      CMAKE_GENERATOR ${gen}
      CMAKE_ARG
        ${LOCAL_CMAKE_BUILD_OPTIONS}
        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
        -DBUILD_TESTING:BOOL=OFF
        -DBUILD_EXAMPLES:BOOL=OFF
        -DITK_LEGACY_REMOVE:BOOL=OFF
        -DModule_MGHIO:BOOL=ON
        -DITKV3_COMPATIBILITY:BOOL=ON #Necessary for IntensityRescaler (niral_utilities)
        -DModule_ITKReview:BOOL=ON #Necessary for MultiAtlasSeg (niral_utilities)
#        -DITK_USE_REVIEW_STATISTICS:BOOL=ON
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/${proj}-install
    )
    FORCE_BUILD_CHECK(${proj})
    set(ITK_DIR ${CMAKE_CURRENT_BINARY_DIR}/${proj}-install/lib/cmake/ITK-4.8)
    set(ITK_DEPEND ${proj}) ## Set the internal dependancy for ITK
  endif()
endmacro()

#-----------------------------------------------------------------------------
# Get and build VTK
#
macro(PACKAGE_NEEDS_VTKWITHQT LOCAL_CMAKE_BUILD_OPTIONS gen)
  set(packageToCheck VTK)
  OPTION(OPT_USE_SYSTEM_${packageToCheck} "Use the system's ${packageToCheck} library." OFF)
  #  MARK_AS_ADVANCED(OPT_USE_SYSTEM_${packageToCheck})
  if(OPT_USE_SYSTEM_VTK)
    find_package(VTK REQUIRED)
    include(${VTK_USE_FILE})
    if( VTK_MAJOR_VERSION VERSION_LESS 5 OR ( VTK_MAJOR_VERSION VERSION_EQUAL 5 AND VTK_MINOR_VERSION VERSION_LESS 6) ) #require vtk5.6 or higher
      message(FATAL_ERROR "VTK 5.6 or more recent is required. Found VTK ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
    endif()
    set(VTK_DEPEND "") ## Set the external depandancy for VTK
  else()
    set(proj vtk-5-6)
    set(vtk_tag -r VTK-5-6)
    set(vtk_module VTK)

    set(vtk_WRAP_TCL OFF)
    set(vtk_WRAP_PYTHON OFF)

    find_package(Qt4 REQUIRED)
    if(QT_USE_FILE)
      include(${QT_USE_FILE})
    endif(QT_USE_FILE)
    set(QT_ARGS
        -DDESIRED_QT_VERSION:STRING=4
        -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
      )

    set(vtk_QT_ARGS)
    set(vtk_QT_ARGS
        ${QT_ARGS}
        -DVTK_USE_GUISUPPORT:BOOL=ON
        -DVTK_USE_QVTK:BOOL=ON
        -DVTK_USE_QT:BOOL=ON
    )
    if(APPLE)
      # Qt 4.6 binary libs are built with empty OBJCXX_FLAGS for mac Cocoa
      set(vtk_QT_ARGS
        ${vtk_QT_ARGS}
        -DVTK_USE_CARBON:BOOL=OFF
        -DVTK_USE_COCOA:BOOL=ON
        -DVTK_REQUIRED_OBJCXX_FLAGS:STRING=
        )
    endif(APPLE)
    set(${proj}_GIT_TAG "v5.6.0")
    set(${proj}_REPOSITORY ${git_protocol}://vtk.org/VTK.git)
    ExternalProject_Add(${proj}
      GIT_REPOSITORY ${${proj}_REPOSITORY}
      GIT_TAG ${${proj}_GIT_TAG}
      UPDATE_COMMAND ""
      SOURCE_DIR ${proj}
      BINARY_DIR ${proj}-build
      CMAKE_GENERATOR ${gen}
      CMAKE_ARGS
        ${LOCAL_CMAKE_BUILD_OPTIONS}
        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
        -DVTK_USE_PARALLEL:BOOL=ON
        -DVTK_DEBUG_LEAKS:BOOL=OFF
        -DVTK_WRAP_TCL:BOOL=${vtk_WRAP_TCL}
        -DVTK_WRAP_PYTHON:BOOL=${vtk_WRAP_PYTHON}
        ${vtk_QT_ARGS}
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/${proj}-install
    )
    FORCE_BUILD_CHECK(${proj})

    set(VTK_DIR ${CMAKE_CURRENT_BINARY_DIR}/${proj}-install/lib/vtk-5.6)
    set(VTK_DEPEND ${proj})
    MESSAGE(STATUS "Setting VTK_DIR to -DVTK_DIR:PATH=${VTK_DIR}")
    set(VTK_CMAKE
       -DVTK_DIR:PATH=${VTK_DIR}
        ${QT_ARGS}
    )
  endif()
endmacro()

macro(PACKAGE_NEEDS_VTK_NOGUI LOCAL_CMAKE_BUILD_OPTIONS gen)
  set(packageToCheck VTK)
  OPTION(OPT_USE_SYSTEM_${packageToCheck} "Use the system's ${packageToCheck} library." OFF)
  #  MARK_AS_ADVANCED(OPT_USE_SYSTEM_${packageToCheck})
  if(OPT_USE_SYSTEM_VTK)
    find_package(VTK 5.6 REQUIRED)
    include(${VTK_USE_FILE})
    set(VTK_DEPEND "") ## Set the external depandancy for ITK
  else()
    set(proj vtk-5-6)

    set(vtk_WRAP_TCL OFF)
    set(vtk_WRAP_PYTHON OFF)

    set(vtk_GUI_ARGS
        -DVTK_USE_GUISUPPORT:BOOL=OFF
        -DVTK_USE_QVTK:BOOL=OFF
        -DVTK_USE_QT:BOOL=OFF
        -DVTK_USE_X:BOOL=OFF
        -DVTK_USE_CARBON:BOOL=OFF
        -DVTK_USE_COCOA:BOOL=OFF
        -DVTK_USE_RENDERING:BOOL=OFF
    )
    if(APPLE)
      # Qt 4.6 binary libs are built with empty OBJCXX_FLAGS for mac Cocoa
      set(vtk_GUI_ARGS
        ${vtk_GUI_ARGS}
        -DVTK_REQUIRED_OBJCXX_FLAGS:STRING=
        )
    endif(APPLE)
    set(${proj}_GIT_TAG "v5.6.0")
    set(${proj}_REPOSITORY ${git_protocol}://vtk.org/VTK.git)
    ExternalProject_Add(${proj}
      GIT_REPOSITORY ${${proj}_REPOSITORY}
      GIT_TAG ${${proj}_GIT_TAG}
      SOURCE_DIR ${proj}
      BINARY_DIR ${proj}-build
      CMAKE_GENERATOR ${gen}
      CMAKE_ARGS
        ${LOCAL_CMAKE_BUILD_OPTIONS}
        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
        -DVTK_USE_PARALLEL:BOOL=ON
        -DVTK_DEBUG_LEAKS:BOOL=OFF
        -DVTK_WRAP_TCL:BOOL=${vtk_WRAP_TCL}
        -DVTK_WRAP_PYTHON:BOOL=${vtk_WRAP_PYTHON}
        ${vtk_GUI_ARGS}
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/${proj}-install
    )
    FORCE_BUILD_CHECK(${proj})

    set(VTK_DIR ${CMAKE_CURRENT_BINARY_DIR}/${proj}-install/lib/vtk-5.6)
    set(VTK_DEPEND ${proj})
    MESSAGE(STATUS "Setting VTK_DIR to -DVTK_DIR:PATH=${VTK_DIR}")
    set(VTK_CMAKE
       -DVTK_DIR:PATH=${VTK_DIR}
    )
  endif()
endmacro()

#-----------------------------------------------------------------------------
# Get and build SlicerExecutionModel
##  Build the SlicerExecutionModel Once, and let all derived project use the same version
macro(PACKAGE_NEEDS_SlicerExecutionModel LOCAL_CMAKE_BUILD_OPTIONS gen)
  set(packageToCheck SlicerExecutionModel)
  OPTION(OPT_USE_SYSTEM_${packageToCheck} "Use the system's ${packageToCheck} library." OFF)
  #  MARK_AS_ADVANCED(OPT_USE_SYSTEM_${packageToCheck})
  if(OPT_USE_SYSTEM_SlicerExecutionModel)
    find_package(GenerateCLP NO_MODULE REQUIRED)
    include(${GenerateCLP_USE_FILE})
    set(GenerateCLP_DEPEND "")
    set(SlicerExecutionModel_DEPEND "")
  else()
    #### ALWAYS BUILD WITH STATIC LIBS
    set(proj SlicerExecutionModel)
    set(${proj}_REPOSITORY "${git_protocol}://github.com/Slicer/SlicerExecutionModel.git")
    set(${proj}_GIT_TAG "e96f2965378cd9a9e5217cd43c8023350138b1cf")
    ExternalProject_Add(${proj}
      GIT_REPOSITORY ${${proj}_REPOSITORY}
      GIT_TAG ${${proj}_GIT_TAG}
      SOURCE_DIR ${proj}
      BINARY_DIR ${proj}-build
      DEPENDS ${ITK_DEPEND}
      CMAKE_GENERATOR ${gen}
      CMAKE_ARGS
        ${LOCAL_CMAKE_BUILD_OPTIONS}
        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
        -DITK_DIR:PATH=${ITK_DIR}
        INSTALL_COMMAND ""
    )
    FORCE_BUILD_CHECK(${proj})
    set(GenerateCLP_DIR ${CMAKE_CURRENT_BINARY_DIR}/${proj}-build/GenerateCLP)
    set(SlicerExecutionModel_DIR ${CMAKE_CURRENT_BINARY_DIR}/${proj}-build)
    set(GenerateCLP_DEPEND "${proj}")
    set(SlicerExecutionModel_DEPEND "${proj}")
  endif()
endmacro()


