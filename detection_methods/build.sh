#!/bin/bash

cd vole
docker build -t phillipecardenuto/cmfd:vole .

cd ../patch_match/
docker build -t phillipecardenuto/cmfd:fd .
