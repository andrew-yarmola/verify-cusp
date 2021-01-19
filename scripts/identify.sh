#!/usr/bin/env bash

cd ../bin

# Test that underflow is detecatble on the system
./test_float

# Run identify
./rootcat data/identify | ./identify

if [ $? -eq 0 ]; then
    echo "Identification succeded."
    # list unique intersecions for ease of checking
    echo "Rerunning to list unique variety intersections."
    ./rootcat data/identify | ./identify | grep "Valid" | sort -u
else
    echo "Identification failed."
fi
