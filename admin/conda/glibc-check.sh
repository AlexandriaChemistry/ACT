#!/bin/bash
# Checks that none of the programs and libraries given as arguments needs a
# newer glibc than that of the sysroot of the conda compiler, which is the
# glibc that the package requires (__glibc). Used by the build scripts of
# the conda packages on Linux, see recipe.yaml.
set -euo pipefail

if [[ $# -eq 0 ]]; then
    echo "No files to check the glibc version of"
    exit 1
fi

SYSROOT=$(${CXX} -print-sysroot)
# features.h separates the words with tabs
LIMIT=$(awk '$1 == "#define" && $2 == "__GLIBC__" {a=$3}
             $1 == "#define" && $2 == "__GLIBC_MINOR__" {b=$3}
             END {if (a != "" && b != "") print a "." b}' \
        "${SYSROOT}/usr/include/features.h")
if [[ -z "${LIMIT}" ]]; then
    echo "No glibc version in ${SYSROOT}/usr/include/features.h"
    exit 1
fi
NEEDED=$(for f in "$@"; do ${OBJDUMP:-objdump} -T "$f"; done |
         grep -oE 'GLIBC_[0-9]+\.[0-9]+' | sed 's/GLIBC_//' | sort -V -u | tail -1)
echo "glibc of the sysroot ${LIMIT}, newest needed ${NEEDED}"
if [[ "$(printf '%s\n%s\n' "${NEEDED}" "${LIMIT}" | sort -V | tail -1)" != "${LIMIT}" ]]; then
    echo "A library needs glibc ${NEEDED}, newer than ${LIMIT} of the sysroot:"
    for f in "$@"; do
        ${OBJDUMP:-objdump} -T "$f" | grep -E "GLIBC_${NEEDED}" | sed "s|^|$f: |"
    done
    exit 1
fi
