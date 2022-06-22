#!/bin/sh

# original author: ramittal@uw.edu, 2020-12-10
# Usage: run_test.sh Invoked from "conda build" to test conda package from meta.yaml
#                    Variables are based on temporary conda env. created and activated
#                    by conda build. All paths used should be in context of conda
#                    special_mat.cpython-37m-x86_64-linux-gnu.so is packaged during
#                    build phase and made available
#                    PREFIX, RECIPE_DIR, SP_DIR, PY_VER are defined by 'conda build'
#                    Not using any custom variables from Makefile, since this is invoked
#                    from conda-build, not from make. conda-build can be invoked by make
# https://docs.conda.io/projects/conda-build/en/latest/user-guide/environment-variables.html
#
# Purpose: Validates presence of limetr files (.so) in python site-packages folder
#          ToDo: Add additional validations to actually invoke limetr code
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

filename="special_mat.cpython-37m-x86_64-linux-gnu.so"

echo "Starting build in working dir :${bname}: at $(date) in '${cwd}'"
echo "Current working directory: '${cwd}'"
echo "Working Conda Environment Variables: PREFIX: '${PREFIX}', RECIPE_DIR='${RECIPE_DIR}'"
# shellcheck disable=SC2154
echo "Package Name: '${PKG_NAME}', Version:'${PKG_VERSION}'"
# shellcheck disable=SC2154
echo "Python version: '${PY_VER}'"
echo "Contents of recipe directory: '${RECIPE_DIR}'"
ls -l "${RECIPE_DIR}/"
echo "Contents of SP_DIR (limetr) dir before test"
ls -l "${SP_DIR}/"

# set expected folder for SO file
SHAREDOBJ_DIR="${SP_DIR}/${PKG_NAME}"
echo "Python's site-packages location: '${SHAREDOBJ_DIR}'"
echo "Target path for package is:'${SHAREDOBJ_DIR}/'"

echo "## Listing conda environment"
conda env list

echo "## Listing installed conda packages"
conda list | tee /tmp/conda_pkg_installed.txt

echo "Validating existence of files at target location"
SHARED_OBJFILE="${SHAREDOBJ_DIR}/${filename}"
echo "Looking for '${SHARED_OBJFILE}'"
if [ -f "${SHARED_OBJFILE}" ] ; then
	echo "File '${SHARED_OBJFILE}' exists at site-packages"
	ls -l "${SHARED_OBJFILE}"
else
	echo "***Error '${SHARED_OBJFILE}' does not exist " >&2
	echo "########## NOT EXITING exit 1"
fi

# validate installation of package
pkgver_file="/tmp/conda-${PKG_NAME}.txt"
grep -i "^${PKG_NAME}" /tmp/conda_pkg_installed.txt > "${pkgver_file}"
exit_status=$?
if [ "0" = "${exit_status}" ] ; then
	# validate version
	version="$(cut -f 2- -d' ' "${pkgver_file}" | xargs | cut -f1 -d' ' | xargs )"
	if [ "${version}" = "${PKG_VERSION}" ] ; then
		echo "Correct version of Package Name: '${PKG_NAME}-${PKG_VERSION}' installed in env."
	else
		echo "***Error Incorrect version of '${PKG_NAME}' is installed in build/test conda env. " >&2
		echo "***Required '${PKG_VERSION}', found: '${version}' " >&2
		exit 1
	fi
else
	echo "***Error '${PKG_NAME}' is not installed in current conda env. " >&2
	exit 1
fi

echo "Checking for existence of tests"
if [ -d "${RECIPE_DIR}/../tests" ] ; then
	test_result=$(python "${RECIPE_DIR}/../tests/check_utils.py" "${RECIPE_DIR}/../tests/check_limetr.py")
	echo "Test result is '${test_result}'"
else
	echo "Unable to find test sources, no test executed"
fi

echo "Completing test of :${PKG_NAME}: at $(date) with exit status: :${exit_status}:"
exit "${exit_status}"

