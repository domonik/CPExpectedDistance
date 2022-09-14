import os
import subprocess
import shutil
import tarfile

FILE_PATH = os.path.abspath(__file__)
MAINDIR = os.path.dirname(os.path.dirname(FILE_PATH))
VRNADIR = os.path.join(MAINDIR, "ViennaRNA")


def build_viennarna(prefix):
    """Builds the ViennaRNA package from source

    """


    # unpack additional sources
    vrna = os.path.join(VRNADIR, "ViennaRNA-2.5.1.tar.gz")
    if os.path.exists(vrna):
        tar = tarfile.open(vrna)
        tar.extractall(VRNADIR)
        tar.close()
    vrna_src = os.path.join(VRNADIR, "ViennaRNA-2.5.1")
    if not os.path.exists(vrna_src):
        raise FileNotFoundError("Did not find ViennaRNA source")
    call = ["./configure", f"--prefix={prefix}"] if prefix is not None else ["./configure"]
    p = subprocess.Popen(
        call, cwd=vrna_src, stderr=subprocess.PIPE, stdout=subprocess.PIPE
    )
    stdout, stderr = p.communicate()
    assert "You can run 'make', 'make check', and 'make install' now!" in stdout.decode()
    p = subprocess.Popen(
        ["make", "-j", str(os.cpu_count())], cwd=vrna_src, stderr=subprocess.PIPE, stdout=subprocess.PIPE
    )
    stdout, stderr = p.communicate()
    p = subprocess.Popen(
        ["make", "install"], cwd=vrna_src, stderr=subprocess.PIPE, stdout=subprocess.PIPE
    )
    stdout, stderr = p.communicate()


if __name__ == '__main__':
    path = "/home/rabsch/miniconda/envs/CPExpectedDistance"
    build_viennarna(prefix=path)
