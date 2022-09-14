from setuptools import setup, find_packages, Extension
import os
import sys
import numpy as np
from sysconfig import get_paths
sys.path.append(os.getcwd())
from tools.build_dependencies import build_viennarna


"-I/home/rabsch/.local/include/ViennaRNA"

CONDA_BUILD = False
SKIP_VIENNARNA = False


def check_vrna(prefix):
    if not os.path.exists(os.path.join(prefix, "lib", "libRNA.a")):
        print("RNAlib missing")
        return False
    return True


filtered_args = []
for i, arg in enumerate(sys.argv):
    if arg == "--conda-build":
        CONDA_BUILD = True
    elif arg == "--skip-viennaRNA":
        SKIP_VIENNARNA = True
    else:
        filtered_args.append(arg)

sys.argv = filtered_args
svi = sys.version_info
python_l = f"python{svi[0]}.{svi[1]}"
py_include = get_paths()["include"]
extra_link_args = []
include_dir = []
if not CONDA_BUILD:
    print("Building ViennaRNA from source")
    prefix = os.getenv("CONDA_PREFIX")
    if prefix is None:
        print("No conda environment found. Falling back to default installation.\n"
              "You might need  root privileges for that\n"
              "Consider using a conda environment instead")
        prefix = "/usr/local"
    _include = os.path.join(prefix, "include")
    _lib = os.path.join(prefix, "lib")
    extra_link_args += [f"-L{_lib}", f"-I{_include}"]
    include_dir += [_include] + [np.get_include()]
    if os.path.exists(os.path.join(prefix, "lib", "libRNA.a")):
        print("assuming ViennaRNA is already built")
        SKIP_VIENNARNA = True

    if not SKIP_VIENNARNA:
        build_viennarna(prefix)
    else:
        print("skipping ViennaRNA build")
    assert check_vrna(prefix), "ViennaRNA was not successfully installed"




extra_link_args += [f"-I{py_include}"]


expected_distance_module = Extension(
    "CPExpectedDistance.c_expected_distance",
    sources=["CPExpectedDistance/c_expected_distance.c"],
    extra_link_args=extra_link_args + ["-lRNA", "-lpthread", "-lstdc++", "-fopenmp", "-lm", f"-l{python_l}", "-Wl,--no-undefined"],
    include_dirs=include_dir,
)



def main():
    setup(
        name="CPExpectedDistance",
        packages=find_packages(),
        version="0.1.0",
        description="Python interface for the fputs C library function",
        author="<your name>",
        author_email="your_email@gmail.com",
        include_package_data=True,
        ext_modules=[expected_distance_module],
        install_requires=[
            "numpy",
        ]
          )

if __name__ == "__main__":
    main()
