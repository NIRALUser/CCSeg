Project(CCSegTool)
set(LOCAL_PROJECT_NAME CCSegTool)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE(ITK REQUIRED)
IF (ITK_FOUND)
   INCLUDE(${USE_ITK_FILE})
ELSE(ITK_FOUND)
   MESSAGE(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
ENDIF (ITK_FOUND)

FIND_PACKAGE(VTK REQUIRED)
IF (VTK_FOUND)
  INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND)
   MESSAGE(FATAL_ERROR, "VTK not found. Please set VTK_DIR.")
ENDIF (VTK_FOUND)

FIND_PACKAGE(GenerateCLP REQUIRED)
IF (GenerateCLP_FOUND)
  INCLUDE(${GenerateCLP_USE_FILE})
ELSE(GenerateCLP_FOUND)
   MESSAGE(FATAL_ERROR, "GenerateCLP not found. Please set GenerateCLP_DIR.")
ENDIF (GenerateCLP_FOUND)

FIND_PACKAGE(Qt4 REQUIRED)
IF(QT_USE_FILE)
  INCLUDE_DIRECTORIES(${QT_INCLUDE_DIR})
  INCLUDE(${QT_USE_FILE})
ELSE(QT_USE_FILE)
   MESSAGE(FATAL_ERROR, "QT not found. Please set QT_DIR.")
ENDIF(QT_USE_FILE)

INCLUDE_DIRECTORIES(
  ${QT_INCLUDE_DIR}
  ${CCSeg_BINARY_DIR}
  ${CCSeg_SOURCE_DIR}
)

QT4_WRAP_CPP(MOC_FILES CCsegtoolguiwindow.h)
QT4_WRAP_UI(UI_FILES CCsegtoolguiwindow.ui)

SET(CCSeg_source CCsegtool.cxx CCsegtoolguiwindow.cxx Imageview.cxx Globalfunc.cxx CCsegtool_computation.cxx CCsegtool_parameters.cxx CCsegtool_initialization.cxx fourier_funcs.cxx VesselRemover.cxx CCcurveProba.cxx ${MOC_FILES} ${UI_FILES})

GENERATECLP(CCsegtool.cxx CCsegtool.xml)
add_executable( CCSeg ${CCSeg_source})

TARGET_LINK_LIBRARIES(CCSeg
  ${ITK_LIBRARIES}
  vtkHybrid
  ${QT_LIBRARIES}
)
