#!/usr/bin/env bash
#Selective initialization of submodules
#Written by by Marko Kosunen, marko.kosunen@aalto.fi, 2017
DIR=$( cd `dirname $0` && pwd )
SUBMODULES="\
    ./thesdk_helpers \
    ./Entities/rtl \
    ./Entities/spice \
    ./Entities/thesdk \
    ./Entities/inverter \
    ./Entities/inv_sim \
    ./Entities/myentity \
    ./Simulations/Slidetemplate \
    ./Simulations/Simtemplate \
    ./Simulations/Inverter"

git submodule sync
for mod in $SUBMODULES; do 
    git submodule update --init $mod
    cd ${mod}
    if [ -f ./init_submodules.sh ]; then
        ./init_submodules.sh
    fi
    cd ${DIR}

done
exit 0


