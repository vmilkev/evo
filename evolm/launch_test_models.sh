#!/bin/bash

uname=$(uname);

case "$uname" in
    Linux)
        echo -n "Linux"
        ;;
    Darwin)
        # adding path to the C++ libplinkio library
        export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/au383883/Documents/MY/codebase/evo/libplinkio/release_lib
        # adding path to the python evolm library/package
        export PYTHONPATH=$PYTHONPATH:/Users/au383883/Documents/MY/codebase/evo/evolm/release
        # adding path to the python evoped library/package
        export PYTHONPATH=$PYTHONPATH:/Users/au383883/Documents/MY/codebase/evo/evoped/release
        #echo "$DYLD_LIBRARY_PATH"
        ;;
    *)
        echo -n "Unknown OS"
        ;;
esac;

# calling the python script
python3 test_models.py