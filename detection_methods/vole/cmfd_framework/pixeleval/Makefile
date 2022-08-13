CC=g++
CFLAGS+=-g
# if your have trouble working with pkg-config, you can also enter the include
# and linking parameters for OpenCV manually. This should read (replace
# /dir/to/opencv/ with your local path) like
#CFLAGS+=-I/dir/to/opencv/include/opencv -I/dir/to/opencv/include
#LDFLAGS+=/dir/to/opencv/lib/libopencv_calib3d.so /dir/to/opencv/lib/libopencv_contrib.so /dir/to/opencv/build/lib/libopencv_core.so
# ... and so on (to be on the safe side, pass all opencv libraries to the
# linker, with their respective path).

CFLAGS+=`pkg-config --cflags opencv`
LDFLAGS+=`pkg-config --libs opencv`

OBJS=pixeleval.o

PROG=pixeleval

$(PROG): $(OBJS)
	$(CC) -o $(PROG) $(OBJS) $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

all: $(PROG)

default: $(PROG)

clean:
	rm -f $(OBJS) $(PROG)

