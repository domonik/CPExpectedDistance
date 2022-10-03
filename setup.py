from setuptools import setup, find_packages, Extension
import os
import sys
import numpy as np
from sysconfig import get_paths
import versioneer
sys.path.append(os.getcwd())
from tools.build_dependencies import build_viennarna



prefix = os.getenv("CONDA_PREFIX")
if prefix is None:
    print("No conda environment found. Falling back to default installation path of ViennaRNA.\n"
          "Consider using a conda environment instead")
    prefix = "/usr/local"


extra_link_args = []
include_dir = []
_include = os.path.join(prefix, "include")
_lib = os.path.join(prefix, "lib")
extra_link_args += [f"-L{_lib}", f"-I{_include}"]
include_dir += [_include] + [np.get_include()]
svi = sys.version_info
python_l = f"python{svi[0]}.{svi[1]}"

if not os.path.exists(os.path.join(prefix, "lib", "libRNA.a")):
    raise FileNotFoundError("Not able to find ViennaRNA RNAlib installation. This version of RNAdist requires ViennaRNA "
                            "to be installed. You can easily install it using Conda:\n"
                            "conda install -c bioconda viennarna"
                            )
if not os.path.exists(os.path.join(prefix, "lib", python_l, "site-packages", "RNA")):
    raise ImportError("Not able to find ViennaRNA python package in your current environment."
                      "Please install it e.g. via Conda\n"
                        "conda install -c bioconda viennarna"
)


svi = sys.version_info
python_l = f"python{svi[0]}.{svi[1]}"
py_include = get_paths()["include"]

extra_link_args += [f"-I{py_include}"]


cp_exp_dist_extension = Extension(
    "CPExpectedDistance.c_expected_distance",
    sources=["CPExpectedDistance/c_expected_distance.c"],
    extra_link_args=extra_link_args + ["-lRNA", "-lpthread", "-lstdc++", "-fopenmp", "-lm", f"-l{python_l}", "-Wl,--no-undefined"],
    include_dirs=include_dir,

)

with open("README.md") as handle:
    LONGDESC = handle.read()


def main():
    setup(
        name="CPExpectedDistance",
        cmdclass=versioneer.get_cmdclass(),
        packages=find_packages(),
        version=versioneer.get_version(),
        description="Python interface for the CP Expected Distances function",
        long_description="Python interface for the CP Expected Distances function",
        author="Domonik",
        author_email="dominik.rabsch@gmail.com",
        include_package_data=True,
        ext_modules=[cp_exp_dist_extension],
        install_requires=[
            "numpy",
        ]
          )

if __name__ == '__main__':
    main()