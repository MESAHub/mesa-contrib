#!/bin/bash

export MESA_CONTRIB_DIR=$(pwd)

do_one(){
    # make workdir for testing
    MESA_WORK_DIR=$(mktemp -d)

    (
        # go there
        cd "${MESA_WORK_DIR}" || exit

        # copy in files from stock workdir
        cp -r "${MESA_DIR}"/star/work/* .

        # copy in example test files
        cp "${MESA_CONTRIB_DIR}"/"$1"/star_example/src/* src/

        # try and make it
        if ./mk > mk.out 2> mk.err; then
            if [ -s 'mk.err' ]; then
                echo "[WARN]"
            else
                echo "[OK]"
            fi
        else
            echo "[FAIL]"
        fi
    ) 
    
    # cleanup
    rm -rf "${MESA_WORK_DIR}"
}

echo
for d in hooks/*/ ; do
    echo "checking ${d%/}"
    do_one "$d"
    echo
done
