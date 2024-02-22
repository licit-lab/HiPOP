# HiPOP
High Performance Optimal Path module

## Install Hipop

Hipop is the C++ graph library used by `mnms`. You must have [CMake](https://cmake.org) (version 3.19 required) and `make` installed on your computer.

### C++ only

Inside your conda environment go to the cpp folder, and install the code using cmake:

```shell
cd cpp
mkdir build
cd build
cmake .. -DCMAKE_PREFIX_PATH=<PREFIX>\
         -DCMAKE_INSTALL_PREFIX=<PREFIX>\
         -DCMAKE_BUILD_TYPE=Release
cmake --build . --target install --config Release
```
Where `<PREFIX>` is the path to your prefix.

If you used conda to install the dependencies, replace it by `$CONDA_PREFIX`.

If you used venv to install the dependencies, replace it by the path to your venv.

You can then lauch the tests in the `build` directory:

```shell
ctest --output-on-failure
```


### Python

To install C++ code use the script `install_cpp.py`:

```shell
python python/install_cpp.py   
```

Then install the python lib:
```shell
python -m pip install python/
```
