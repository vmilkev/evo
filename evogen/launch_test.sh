#!/bin/bash

uname=$(uname);

case "$uname" in
    Linux)
        echo -n "Linux"
        # adding evogen library to the python path
        export PYTHONPATH=$PYTHONPATH:/Users/au383883/Documents/MY/codebase/evo/evogen/release
        ;;
    Darwin)
        # adding evogen library to the python path
        export PYTHONPATH=$PYTHONPATH:/Users/au383883/Documents/MY/codebase/evo/evogen/release
        ;;
    *)
        echo -n "Unknown OS"
        ;;
esac;

# calling the python script
python3 tests/data/discrete_gen/anc_pop2_with_comments.py