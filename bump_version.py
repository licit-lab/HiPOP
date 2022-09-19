import re
from pathlib import Path
import argparse

PATTERN_VERSION = r"\d+\.\d+\.\d+"

files = [("cpp/CMakeLists.txt", re.compile(rf"(project\(HiPOP VERSION )({PATTERN_VERSION})(\))")),
         ("python/hipop/__init__.py", re.compile(rf'(__version__ = ")({PATTERN_VERSION})(")')),
         ("conda/recipe/meta.yaml", re.compile(r'({% set version = ")('+PATTERN_VERSION+r')(" %})'))]


def read(file):
    with open(file, "r") as f:
        contents = f.read()
    return contents


def write(file, new_contents):
    with open(file, "w") as f:
        f.write(new_contents)


cwd = Path(__file__).parent.resolve()


def show_diff(contents, new_contents):
    for i, (lc, lnc) in enumerate(zip(contents.split("\n"), new_contents.split("\n"))):
        if lc != lnc:
            print(f"\tline {i}: {lc} -> {lnc}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('MAJOR', help='major version number', type=int)
    parser.add_argument('MINOR', help='minor version number', type=int)
    parser.add_argument('PATCH', help='patch version number', type=int)
    args = parser.parse_args()

    new_version = f"{args.MAJOR}.{args.MINOR}.{args.PATCH}"

    for f, pattern in files:
        print("Bumping version in file:", f)
        f_path = cwd.joinpath(f)
        contents = read(f_path)
        new_contents = pattern.sub(f"\g<1>{new_version}\g<3>", contents, count=1)
        show_diff(contents, new_contents)
        write(f_path, new_contents)
