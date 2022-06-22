# check if the requried the packages are installed
import os
import pathlib
import importlib
from sys import platform


# installation of the required packages
required_modules = [('numpy',
                     'conda install -y numpy'),
                    ('scipy',
                     'conda install -y scipy'),
                    ('ipopt',
                     'conda install -y -c conda-forge cyipopt')]


def check_module(module_name, install_command):
    try:
        importlib.import_module(module_name)
    except:
        os.system(install_command)


def extract_lib(lib_name, des_lib_folder):
    conda_lib = os.path.join(os.getenv("CONDA_PREFIX"), "lib")
    pathlib.Path(des_lib_folder).mkdir(exist_ok=True)
    lib_files = [file_name for file_name in os.listdir(conda_lib)
                 if lib_name in file_name]

    if not lib_files:
        raise FileNotFoundError(lib_name + "not found!")

    for file in lib_files:
        os.system(" ".join(["cp -L",
                            os.path.join(conda_lib, file),
                            des_lib_folder]))

    if platform == "linux" or platform == "linux2":
        required_lib_name = lib_name + ".so"
        related_lib_files = [file_name for file_name in lib_files
                             if required_lib_name in file_name]
        assert any(related_lib_files)
        if not pathlib.Path(os.path.join(des_lib_folder,
                                         required_lib_name)).exists():
            os.system(" ".join(["ln -s",
                                related_lib_files[-1],
                                os.path.join(des_lib_folder,
                                             required_lib_name)]))


for module_name, install_command in required_modules:
    check_module(module_name, install_command)

# create the library of blas and lapack
des_lib_folder = "./lib"
extract_lib("libblas", des_lib_folder)
extract_lib("liblapack", des_lib_folder)
