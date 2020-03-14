#!/bin/bash

export MESA_CONTRIB_DIR=$(pwd)

do_one(){
    # make workdir for testing
    MESA_WORK_DIR=$(mktemp -d)

    (
        # go there
        cd "${MESA_WORK_DIR}"

        # copy in files from stock workdir
        cp -r "${MESA_DIR}"/star/work/* .

        # copy in example test files
        cp "${MESA_CONTRIB_DIR}"/"$1"/star_example/src/* src/

        # try and make it
        if ./mk > mk.out; then
            echo "  [OK]"
        else
            echo "  [FAIL]"
        fi
    ) 
    
    # cleanup
    rm -rf "${MESA_WORK_DIR}"
}

for d in hooks/*/ ; do
    echo
    echo "****"
    echo "checking $d"
    do_one "$d"
    echo
done
