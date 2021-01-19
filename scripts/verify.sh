#!/usr/bin/env bash

cd ../bin

# Test that underflow is detecatble on the system
./test_float

# Run identify
./rootcat data/verify | ./verify

if [ $? -eq 0 ]; then
    echo "Verification succeded."
else
    echo "Verification failed."
fi
