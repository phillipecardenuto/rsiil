add_library(utilities
	Utilities/iio.c
	Utilities/LibImages.cpp
	Utilities/Utilities.cpp
	Utilities/LibMatrix.cpp
	Utilities/filters.cpp
	Utilities/FeatManager/zMManager.cpp
	Utilities/FeatManager/siftManager.cpp
	Utilities/PatchMatch/patchmatch.cpp
	)

add_library(vlfeat
	vlfeat/vl/host.c
	vlfeat/vl/generic.c
	vlfeat/vl/dsift.c
	vlfeat/vl/scalespace.c
	vlfeat/vl/sift.c
	vlfeat/vl/random.c
	vlfeat/vl/imopv.c
	vlfeat/vl/imopv_sse2.c
	)

include_directories(vlfeat)

add_executable(forgery_detection main.cpp)
target_link_libraries(forgery_detection utilities vlfeat ${LAPACK_LIBRARIES} ${CBLAS_LIBRARIES})
