#!/bin/sh

# original author: ramittal@uw.edu, 2020-12-10
# Usage: build.sh    Invoked from "conda build" to build conda package from meta.yaml
#                    ../build/lib.linux-x86_64-3.7/limetr/special_mat.cpython-37m-x86_64-linux-gnu.so
#                    must have been build and available before invoking this script
#                    PREFIX, RECIPE_DIR, SP_DIR, PY_VER are defined by 'conda build'
#                    Not using any custom variables from Makefile, since this is invoked
#                    from conda-build, not from make. conda-build can be invoked by make
# https://docs.conda.io/projects/conda-build/en/latest/user-guide/environment-variables.html
#
# Purpose: Produces conda package including .so file for limetr
#

cwd="$(pwd)"
bname="$(basename "${cwd}" )"

# check and display environment variables used, never null from conda-build
# shellcheck disable=SC2154
if [ "" = "${RECIPE_DIR}" ] ; then
	echo "RECIPE_DIR is required but found empty"
	exit 1
else
	echo "RECIPE_DIR is '${RECIPE_DIR}'"
fi
# shellcheck disable=SC2154
if [ "" = "${PREFIX}" ] ; then
	echo "PREFIX is required but found empty"
	exit 1
else
	echo "PREFIX is '${PREFIX}'"
fi
# shellcheck disable=SC2154
if [ "" = "${SP_DIR}" ] ; then
	echo "SP_DIR is required but found empty"
	exit 1
else
	echo "SP_DIR is '${SP_DIR}'"
fi

SOURCE_DIR="${RECIPE_DIR}/../src/limetr"
filename="special_mat.cpython-37m-x86_64-linux-gnu.so"
SHARED_OBJFILE="${RECIPE_DIR}/../build/lib.linux-x86_64-3.7/limetr/${filename}"

# set target folder for SO file
# shellcheck disable=SC2154
TARGET_DIR="${SP_DIR}/${PKG_NAME}"

echo "Starting build in working folder :${bname}: at $(date) in '${cwd}'"
echo "Current working directory: '${cwd}'"
echo "Working Conda Environment Variables: PREFIX: '${PREFIX}', RECIPE_DIR='${RECIPE_DIR}'"
echo "Python's site-packages location: '${TARGET_DIR}'"
# shellcheck disable=SC2154
echo "Package Name: '${PKG_NAME}', Version:'${PKG_VERSION}'"
# shellcheck disable=SC2154
echo "Python version: '${PY_VER}'"
echo "Contents of recipe directory: '${RECIPE_DIR}'"
ls -l "${RECIPE_DIR}/"
echo "Target path for package is:'${TARGET_DIR}/'"
echo "Contents of SP_DIR (limetr) dir before copy"
ls -l "${SP_DIR}/"

# validate valid package and version names
if [ "none" = "${PKG_NAME}" ] ||  [ "None" = "${PKG_VERSION}" ] ; then
	echo "***Error Invalid Package or version values " >&2
	echo 1
fi

echo "Starting copy of artifacts to '${TARGET_DIR}/"
if [ -d "${TARGET_DIR}" ] ; then
	echo "Target dir '${TARGET_DIR} exists"
else
	echo "Creating not existing target dir '${TARGET_DIR}"
	mkdir -p "${TARGET_DIR}"
fi
# copy shared library to target location
if [ -f "${SHARED_OBJFILE}" ] ; then
	echo "Copying '${SHARED_OBJFILE}'"
	cp -f "${SHARED_OBJFILE}" "${TARGET_DIR}/"
	exit_status="$?"
	if [ "${exit_status}" -ne 0 ] ; then
		echo "***Error '${exit_status}'*** during copy of '${SHARED_OBJFILE}' file into site-packages dir " >&2
		exit 1
	fi
else
	echo "***Error '${SHARED_OBJFILE}' does not exist at source " >&2
	exit 1
fi
echo "Copying source files from '${SOURCE_DIR}' to '${TARGET_DIR}'"
ls -l "${SOURCE_DIR}"
for afile in "${SOURCE_DIR}"/*.py
do
	echo "Copying ${afile}"
	cp -f "${afile}" "${TARGET_DIR}/".
	exit_status="$?"
	if [ "0" = "${exit_status}" ] ; then
		echo "${afile} copied successfully"
	else
		echo "*** Copy of ${afile} failed with status: '${exit_status}' into '${TARGET_DIR}' dir " >&2
		exit 1
	fi
done

echo "Contents of SP_DIR (limetr) dir after copy"
ls -l "${SP_DIR}"
if [ -d "${TARGET_DIR}" ] ; then
	echo "Contents of '${TARGET_DIR}'"
	ls -l "${TARGET_DIR}"
else
	echo "*** Error: expected '${TARGET_DIR}' does not exist"
	exit 1
fi

echo "Completing build of :${PKG_NAME}: at $(date)"
exit "0"

