#!/usr/bin/env bash

# build
cd ../src
make test_float
make rootcat
make verify
make identify
