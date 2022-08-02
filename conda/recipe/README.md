# Build conda pkgs

You need to have the conda package `conda-build` installed in your `base` environment.

At the root of `hipop` run the following command inside a terminal:

```shell
conda build . -c conda-forge --output-folder conda-bld
```

If it is a success, you will find the package inside the `conda-bld` folder.