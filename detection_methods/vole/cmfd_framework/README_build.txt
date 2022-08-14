To compile the code, you require
- OpenCV (tested with 2.4.0)
- Boost (tested with 1.45)
- cmake (tested with 2.8.2)

We typically built the code on a debian squeeze Linux. The build process worked
also without modifications under Ubuntu Linux.

enter the root directory of the code (cmfd_framework). To build the code, execute these steps:
----------------snip------------------------
mkdir build
cd build
ccmake ../
----------------snap------------------------

The curses interface of cmake shows up. Set the boolean switches
Vole_CMFD = ON
Vole_CMFD_Ground_Truth = ON
Vole_Shell = ON

Set the variable OpenCV_DIR to the directory in your OpenCV installation that
contains the file OpenCVConfig.cmake . This is typically <install_dir>/share/OpenCV/

If you want to use a custom build of boost, press 't' (for detailed options) and fix the paths to the boost libraries, namely
Boost_DATE_TIME_LIBRARY
...
Boost_THREAD_LIBRARY_RELEASE
and also
Boost_INCLUDE_DIR
Boost_LIBRARY_DIRS
If your system provides a reasonably recent version of boost, this should not be necessary.

Press 'c' to configure the build.

If everything went well, press 'g' to create Makefiles and to leave ccmake.

Type 'make' to build the code.

If everything goes well, you should have a binary bin/vole.
You will find a description on how to use the binary in README_commands.txt


