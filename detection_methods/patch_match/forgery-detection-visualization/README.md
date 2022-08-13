% Copy move forgery detection.

# ABOUT

* Author    : EHRET Thibaud <ehret.thibaud@gmail.com>
* Copyright : (C) 2017 IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see GPLv3.txt

Version 1.0, released on June 19, 2017

# OVERVIEW

This source code provides an implementation of a copy move forgery detector originally developped in "Davide Cozzolino, Giovanni Poggi, and Luisa Verdoliva, Efficient dense-field copy-move forgery detection, IEEE Transactions on Information Forensics and Security, 10 (2015),pp. 2284-2297", which has been studied in details in http//www.ipol.im/ (TODO give the exact link to the article).

# UNIX/LINUX/MAC USER GUIDE

The code is compilable on Unix/Linux and Mac OS. 

- Compilation. 
Automated compilation requires the cmake and make programs.

- Library. 
This code requires the libpng, libtiff, libjpeg and cblas library (libopenblas-dev).

- Image format. 
Only the TIFF, JPEG and PNG formats are supported. 
 
-------------------------------------------------------------------------
Usage:
1. Download the code package and extract it. Go to that directory. 

2. Configure and compile the source code using cmake and make. 
It is recommended that you create a folder for building:

UNIX/LINUX:
$ mkdir build; cd build
$ cmake ..
$ make

MAC/Homebrew:
$ mkdir build; cd build
$ brew install openblas
$ cmake .. -DCBLAS_INCLUDES=$(brew --prefix openblas)/include -DCBLAS_LIBRARIES=$(brew --prefix openblas)/lib/libblas.dylib
$ make

Binaries will be created in build/bin folder.

NOTE: By default, the code is compiled with OpenMP multithreaded
parallelization enabled (if your system supports it). 
The number of threads used by the code is defined in Utilities/parameters.h.

3. Usage instruction:
Running "./forgery_detection -i forged.png" computes the result using the method and parameters from "Davide Cozzolino, Giovanni Poggi, and Luisa Verdoliva, Efficient dense-field copy-move forgery detection, IEEE Transactions on Information Forensics and Security, 10 (2015),pp. 2284-2297".
"./forgery_detection --help" list all available options

4. A comprehensive result is available in "result.txt".
The algorithm also produces different images presenting the result after a specific test:
"filteredPmDisp.png" corresponds to a visual representation of the displacement produced by PatchMatch after median filtering
"errorMap.png" corresponds to a visual representation of the error computed based on the the displacement map. The whiter the image, the higher the error
"detectionMask.png" corresponds to the initial mask after simple thresholding
"filteredMask.png" corresponds to final mask after all post-processing

5. An example of forged image used in the example section of the joint article (rotation) is available in "example" as "forged.png". A possible result of the run of this program is also available in the same folder. The true forgery is shown in "original.png". This example has been generated using the command './forgery_detection -i forged.png'.

6. This project contains the following source files on top of the sources of VLFeat (http://www.vlfeat.org/):
	src/CMakeLists.txt
	src/main.cpp
	src/Utilities/iio.h
	src/Utilities/filters.cpp
	src/Utilities/filters.h
	src/Utilities/parameters.h
	src/Utilities/iio.c
	src/Utilities/LibImages.h
	src/Utilities/PatchMatch/patchmatch.cpp
	src/Utilities/PatchMatch/patchmatch.h
	src/Utilities/cmd_option.h
	src/Utilities/Utilities.cpp
	src/Utilities/LibImages.cpp
	src/Utilities/Utilities.h
	src/Utilities/LibMatrix.h
	src/Utilities/tools.h
	src/Utilities/FeatManager/zMManager.h
	src/Utilities/FeatManager/featManager.h
	src/Utilities/FeatManager/siftManager.cpp
	src/Utilities/FeatManager/zMManager.cpp
	src/Utilities/FeatManager/siftManager.h
	src/Utilities/LibMatrix.cpp


7. The files that have been reviewed for IPOL publication are
	src/main.cpp
	src/Utilities/filters.h
	src/Utilities/filters.cpp
	src/Utilities/PatchMatch/patchmatch.cpp
	src/Utilities/PatchMatch/patchmatch.h
	src/Utilities/FeatManager/zMManager.h
	src/Utilities/FeatManager/zMManager.cpp


# ABOUT THIS FILE

Copyright 2017 IPOL Image Processing On Line http://www.ipol.im/

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.
