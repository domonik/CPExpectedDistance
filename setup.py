from distutils.core import setup, Extension


"-I/home/rabsch/.local/include/ViennaRNA"


module_1 = Extension(
    "exp_d",
    sources=["CPExpectedDistance/c_expected_distance.c"],
    extra_link_args=["-lRNA", "-lpthread", "-lstdc++", "-fopenmp", "-lm", "-Wl,--no-undefined"],
    # include_dirs=["/home/rabsch/.local/include/", "/home/rabsch/.local/lib/"],
)


def main():
    setup(
        name="cp_expected_distance",
        version="0.1.0",
        description="Python interface for the fputs C library function",
        author="<your name>",
        author_email="your_email@gmail.com",
        include_package_data=True,

        ext_modules=[module_1]

          )

if __name__ == "__main__":
    main()
