import os
import subprocess
import sys
import tempfile
from pathlib import Path


def build_hipop_cpp(prefix):
    src_dir = Path(__file__).parent.parent.joinpath("cpp").resolve()
    temp_dir = tempfile.TemporaryDirectory()
    os.chdir(temp_dir.name)
    try:
        subprocess.run([fr"cmake -S {src_dir} -DCMAKE_PREFIX_PATH={prefix} -DCMAKE_INSTALL_PREFIX={prefix} -DCMAKE_BUILD_TYPE=Release"], shell=True)
        subprocess.run([r"cmake --build . --target install --config Release"], shell=True)
    except:
        print("Build failed ..")
        sys.exit(-1)
    finally:
        temp_dir.cleanup()
        os.chdir(Path(__file__).parent)


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




