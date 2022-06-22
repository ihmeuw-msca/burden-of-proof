# make file for pynlme class
OUTPUT_DIR=output
CONDA_PKG_DIR=conda_pkg
NUMPY_VER=1.19.1

.PHONY: clean, tests

build: setup.py
	python setup.py build

install: setup.py
	python check_requirements.py
	python setup.py install

sdist: setup.py
	python setup.py sdist

tests:
	python tests/check_utils.py
	python tests/check_limetr.py

clean:
	find . -name "*.so*" | xargs rm -rf
	find . -name "*.pyc" | xargs rm -rf
	find . -name "__pycache__" | xargs rm -rf
	find . -name "build" | xargs rm -rf
	find . -name "dist" | xargs rm -rf
	find . -name "MANIFEST" | xargs rm -rf
	rm -rf ./lib ./$(OUTPUT_DIR)

uninstall:
	find $(CONDA_PREFIX)/lib/ -name "*limetr*" | xargs rm -rf

package: src/limetr/Makefile src/limetr/special_mat.f90 src/limetr/utils.py
	@echo "### Ensure version number in $(CONDA_PKG_DIR) matches with version in setup.py"
	@echo "Currently only tested for linux"
	@echo "Installing conda pre-requirements"
	@conda install --yes --strict-channel-priority -c conda-forge -c defaults conda-build conda-verify
	@echo "Installing additional conda dependencies"
	@conda install --yes --strict-channel-priority -c conda-forge -c defaults numpy==1.19.1 scipy==1.5.2 cyipopt
	@echo "Building conda package for limetr (from $(CONDA_PKG_DIR) folder)"
	conda build -k --no-anaconda-upload --verify --numpy $(NUMPY_VER) --output-folder "$(OUTPUT_DIR)" --cache-dir /tmp/limetrcache ./$(CONDA_PKG_DIR)/
	@echo "conda build status:'$?'"
	@echo "Generated conda package file is: $(OUTPUT_DIR)/linux-64/limetr*.tar.bz2"
	@ls -l $(OUTPUT_DIR)/linux-64/limetr*.tar.bz2

