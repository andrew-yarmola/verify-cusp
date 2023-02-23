#!/usr/bin/env bash

pushd ../src
mkdir -p ../bin
# build
make tests
make rootcat
make verify
make identify
popd
