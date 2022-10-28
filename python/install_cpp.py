import os
import subprocess
import sys
import tempfile
from pathlib import Path
import platform

def build_hipop_cpp(prefix):
    cwd = os.getcwd()
    src_dir = Path(__file__).parent.parent.joinpath("cpp").resolve()
    temp_dir = tempfile.TemporaryDirectory()
    os.chdir(temp_dir.name)
    try:
        subprocess.check_call(fr"cmake -S '{src_dir}' -DCMAKE_PREFIX_PATH={prefix} -DCMAKE_INSTALL_PREFIX={prefix} -DCMAKE_BUILD_TYPE=Release", shell=True)
        subprocess.check_call(r"cmake --build . --target install --config Release", shell=True)
    except:
        print("Build failed ..")
        sys.exit(-1)
    finally:
        if platform.system() != "Windows":
            temp_dir.cleanup()
        os.chdir(cwd)


def main():
    try:
        subprocess.check_output(['cmake', '--version'], stderr=subprocess.STDOUT).decode('utf-8')
    except:
        print("Cmake not found, please install it")
        sys.exit(-1)

    print("Compiling and installing C++ sources ...")
    build_hipop_cpp(sys.prefix)
    print(f"Successfully installed in: {sys.prefix}")


if __name__ == "__main__":
    main()




