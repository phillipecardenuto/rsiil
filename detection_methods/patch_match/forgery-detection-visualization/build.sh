mkdir build
cd build
cmake ../
make -j
cp -f bin/forgery_detection ../
cd ../
rm -fr build