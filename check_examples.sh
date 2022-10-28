#!/bin/bash

# export environment variables that MESA needs
export MESA_CONTRIB_DIR=$(pwd)
export MESA_FORCE_PGSTAR_FLAG=false

# make tmpdir for testing
# local to script; MESA doesn't need it so no need to export
MESA_WORK_DIR_BASE=$(mktemp -d)

# arguments are
# 1. folder name of hook
# 2. module name (star/binary/astero)
do_one(){

    (
	echo checking "$1" "$2"

        EXAMPLE_DIR="${MESA_CONTRIB_DIR}"/"$1""$2"_example

        # if this kind of example doesn't exist, exit
        if [ ! -d "${EXAMPLE_DIR}" ]; then
	    exit
	fi

        # set MESA_WORK_DIR
        MESA_WORK_DIR="${MESA_WORK_DIR_BASE}"/"$1""$2"
        mkdir -p "${MESA_WORK_DIR}"

        # go there
        cd "${MESA_WORK_DIR}" || exit

        # copy in files from stock workdir
        cp -r "${MESA_DIR}"/"$2"/work/* .

        # copy in example test files
        # cp -vr "${EXAMPLE_DIR}"/src/* src/
        cp -vr "${EXAMPLE_DIR}"/* ./

        # try and make it
        if ./mk > mk.out 2> mk.err; then
            if [ -s 'mk.err' ]; then
                echo "[COMPILE WARN]"
                echo "see $MESA_WORK_DIR for details"
                exit
            else
                echo "[COMPILE OK]"
                # cleanup
                # rm -rf "${MESA_WORK_DIR}"
            fi
        else
            echo "[COMPILE FAIL]"
            echo "see $MESA_WORK_DIR for details"
            exit
        fi

	# try to run it
	if ./rn > rn.out 2> rn.err; then
            if [ -s 'rn.err' ]; then
		echo "[RUN WARN]"
                echo "see $MESA_WORK_DIR for details"
                exit
            else
                echo "[RUN OK]"
                # cleanup
                rm -rf "${MESA_WORK_DIR}"
            fi
        else
            echo "[RUN FAIL]"
            echo "see $MESA_WORK_DIR for details"
            exit
        fi
    ) 
    
}

echo
for d in hooks/*/ ; do
    do_one "$d" star
    do_one "$d" binary
    do_one "$d" astero
    echo
done
