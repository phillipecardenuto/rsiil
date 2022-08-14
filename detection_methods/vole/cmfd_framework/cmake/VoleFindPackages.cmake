# Library versions
set(VOLE_MINIMUM_OPENCV_VERSION "2.2.0")
set(VOLE_MINIMUM_QT_VERSION "4.6.0")
set(VOLE_MINIMUM_BOOST_VERSION "1.35")
set(VOLE_MINIMUM_EIGEN_VERSION "3.0")
set(VOLE_MINIMUM_SLEPC_VERSION "3.1") # TODO: version check missing in FindSLEPc

# OpenCV
find_package(OpenCV PATHS "/net/cv/lib/share/OpenCV" "/local/opencv/share/OpenCV")
if (${OpenCV_VERSION})
	if(${OpenCV_VERSION} VERSION_LESS ${VOLE_MINIMUM_OPENCV_VERSION})
		set(OpenCV_FOUND)
	endif()
endif()
if(OpenCV_FOUND)
	add_definitions(-DOPENCV_VERSION=${OpenCV_VERSION})
	add_definitions(-DWITH_OPENCV2)
endif()
vole_check_package(OPENCV
	"OpenCV"
	"Please install OpenCV >=${VOLE_MINIMUM_OPENCV_VERSION} or set OpenCV_DIR."
	OpenCV_FOUND
	"${OpenCV_INCLUDE_DIRS}"
	"${OpenCV_LIBS}"
)

# Thread Building Blocks
find_package(TBB)
vole_check_package(TBB
	"TBB"
	"Please install TBB"
	TBB_FOUND
	"${TBB_INCLUDE_DIRS}"
	"${TBB_LIBRARIES}"
)

# OpenGL
find_package(OpenGL)
vole_check_package(OPENGL
	"OpenGL"
	"Please check your OpenGL installation."
	OpenGL_FOUND
	"${OpenGL_INCLUDE_DIR}"
	"${OpenGL_LIBRARIES}"
)

# QtGui
find_package(Qt4 ${VOLE_MINIMUM_QT_VERSION} COMPONENTS QtCore QtGui)
vole_check_package(QT
	"Qt4"
	"Please install Qt4 >=${VOLE_MINIMUM_QT_VERSION} or set QT_QMAKE_EXECUTABLE."
	QT_FOUND
	"${QT_INCLUDE_DIR};${QT_QTCORE_INCLUDE_DIR};${QT_QTGUI_INCLUDE_DIR}"
	"${QT_QTCORE_LIBRARY};${QT_QTGUI_LIBRARY}"
)

# QtOpenGL
find_package(Qt4 ${VOLE_MINIMUM_QT_VERSION} COMPONENTS QtOpenGL)
vole_check_package(QT_OPENGL
	"Qt4 OpenGL"
	"Please install Qt4 >=${VOLE_MINIMUM_QT_VERSION} or set QT_QMAKE_EXECUTABLE."
	QT_QTOPENGL_FOUND
	"${QT_INCLUDE_DIR};${QT_QTOPENGL_INCLUDE_DIR}"
	"${QT_QTOPENGL_LIBRARY}"
)

# QtXml
find_package(Qt4 ${VOLE_MINIMUM_QT_VERSION} COMPONENTS QtXml)
vole_check_package(QT_XML
	"Qt4 XML"
	"Please install Qt4 >=${VOLE_MINIMUM_QT_VERSION} or set QT_QMAKE_EXECUTABLE."
	QT_QTXML_FOUND
	"${QT_INCLUDE_DIR};${QT_QTXML_INCLUDE_DIR}"
	"${QT_QTXML_LIBRARY}"
)

# ITK
find_package(ITK)
vole_check_package(ITK
    "ITK"
    "Please install ITK."
    ITK_FOUND
    "${ITK_INCLUDE_DIRS}"
    "${ITK_LIBRARIES}"
)

# Boost
if(WIN32)
	set(Boost_USE_STATIC_LIBS ON)
	set(BOOST_ROOT "C:\\boost" CACHE STRING "Boost Root Directory.")
else()
	set(BOOST_ROOT "" CACHE STRING "Boost Root Directory.")
endif()
find_package(Boost ${VOLE_MINIMUM_BOOST_VERSION})
# if you need the boost version in the code, #include <boost/version.hpp>
#if(Boost_FOUND)
#	add_definitions(-DBOOST_VERSION=${Boost_VERSION})
#endif()
vole_check_package(BOOST
	"Boost"
	"Please install Boost system >= ${VOLE_MINIMUM_BOOST_VERSION} or set Boost_ROOT."
	Boost_FOUND
	"${Boost_INCLUDE_DIR}/include/;${Boost_INCLUDE_DIR}"
	""
)

# Boost system
find_package(Boost ${VOLE_MINIMUM_BOOST_VERSION} COMPONENTS system filesystem program_options serialization thread date_time python)
vole_check_package(BOOST_SYSTEM
	"Boost system"
	"Please install Boost system >= ${VOLE_MINIMUM_BOOST_VERSION} or set Boost_ROOT."
	Boost_SYSTEM_FOUND
	"${Boost_INCLUDE_DIR}/include/;${Boost_INCLUDE_DIR}"
	"${Boost_SYSTEM_LIBRARY}"
)

# Apply checks for the specified packages
# Boost filesystem
vole_check_package(BOOST_FILESYSTEM
	"Boost filesystem"
	"Please install Boost filesystem >= ${VOLE_MINIMUM_BOOST_VERSION} or set Boost_ROOT."
	Boost_FILESYSTEM_FOUND
	"${Boost_INCLUDE_DIR}/include/;${Boost_INCLUDE_DIR}"
	"${Boost_FILESYSTEM_LIBRARY}"
)
# Boost program options
vole_check_package(BOOST_PROGRAM_OPTIONS
	"Boost program options"
	"Please install Boost program options >=${VOLE_MINIMUM_BOOST_VERSION} or set Boost_ROOT."
	Boost_PROGRAM_OPTIONS_FOUND
	"${Boost_INCLUDE_DIR}/include/;${Boost_INCLUDE_DIR}"
	"${Boost_PROGRAM_OPTIONS_LIBRARY}"
)
# Boost serialization
vole_check_package(BOOST_SERIALIZATION
	"Boost Serialization"
	"Please install Boost serialization >=${VOLE_MINIMUM_BOOST_VERSION} or set Boost_ROOT."
	Boost_SERIALIZATION_FOUND
	"${Boost_INCLUDE_DIR}/include/;${Boost_INCLUDE_DIR}"
	"${Boost_SERIALIZATION_LIBRARY}"
)
# Boost thread
vole_check_package(BOOST_THREAD
	"Boost thread"
	"Please install Boost thread >=${VOLE_MINIMUM_BOOST_VERSION} or set Boost_ROOT."
	Boost_THREAD_FOUND
	"${Boost_INCLUDE_DIR}/include/;${Boost_INCLUDE_DIR}"
	"${Boost_THREAD_LIBRARY}"
)
# Boost date time
vole_check_package(BOOST_DATE_TIME
	"Boost date time"
	"Please install Boost date time >=${VOLE_MINIMUM_BOOST_VERSION} or set Boost_ROOT."
	Boost_DATE_TIME_FOUND
	"${Boost_INCLUDE_DIR}/include/;${Boost_INCLUDE_DIR}"
	"${Boost_DATE_TIME_LIBRARY}"
)
# Boost python
vole_check_package(BOOST_PYTHON
	"Boost python"
	"Please install Boost python >=${VOLE_MINIMUM_BOOST_VERSION} or set Boost_ROOT."
	Boost_PYTHON_FOUND
	"${Boost_INCLUDE_DIR}/include/;${Boost_INCLUDE_DIR}"
	"${Boost_PYTHON_LIBRARY}"
)

# Python
find_package(PythonLibs)
vole_check_package(PYTHON
	"Python"
	"Please install Python"
	PYTHONLIBS_FOUND
	"${PYTHON_INCLUDE_DIR}"
	"${PYTHON_LIBRARIES}"
)

# Numpy
find_package(NumPy)
vole_check_package(NUMPY
	"Numpy"
	"Please install Numpy"
	NUMPY_FOUND
	"${NUMPY_INCLUDE_DIRS}"
	""
)

# PETSc && SLEPc
find_package(PETSc)
vole_check_package(PETSC
	"PETSc"
	"Please install PETSc (and SLEPc) or set PETSC_DIR."
	PETSC_FOUND
	"${PETSC_INCLUDE_DIR}"
	"${PETSC_LIBRARIES}"
)
find_package(SLEPc)
vole_check_package(SLEPC
	"SLEPc"
	"Please install SLEPc >=${VOLE_MINIMUM_SLEPC_VERSION} or set SLEPC_DIR."
	SLEPC_FOUND
	"${SLEPC_INCLUDE_DIR}"
	"${SLEPC_LIBRARIES}"
)

# libEigen3
find_package(Eigen3 ${VOLE_MINIMUM_EIGEN_VERSION})
vole_check_package(EIGEN3
	"Eigen"
	"Please install libeigen3 >=${VOLE_MINIMUM_EIGEN_VERSION}"
	EIGEN3_FOUND
	"${EIGEN3_INCLUDE_DIR}"
	""
)

# flann
find_package(Flann)
vole_check_package(FLANN
        "Flann"
        "Please install flann"
        FLANN_FOUND
        "${FLANN_INCLUDE_DIRS}"
        "${FLANN_LIBRARIES}"
)

