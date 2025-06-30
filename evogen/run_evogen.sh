#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <python_script>"
    exit 1
fi

uname=$(uname);

case "$uname" in
    Linux)
        echo -n "Linux"
        # adding evogen library to the python path
        export PYTHONPATH=$PYTHONPATH:/home/vimi/Documents/codebase/evo/evogen/release
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
python3 $1
