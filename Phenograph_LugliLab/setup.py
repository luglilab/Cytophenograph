from setuptools import setup, find_packages
import sys


if sys.version_info.major != 3:
    raise RuntimeError("PhenoGraph requires Python 3")

main_ns = {}

# get version
with open("phenograph/version.py") as f:
    exec(f.read(), main_ns)
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="PhenoGraph",
    description="Graph-based clustering for high-dimensional single-cell data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=main_ns["__version__"],
    author=main_ns["__author__"],
    author_email=main_ns["__email__"],
    packages=find_packages(),
    package_data={
        "": ["louvain/*convert*", "louvain/*community*", "louvain/*hierarchy*"]
    },
    include_package_data=True,
    zip_safe=False,
    url="https://github.com/dpeerlab/PhenoGraph.git",
    license="LICENSE",
    install_requires=open("requirements.txt").read(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Operating System :: POSIX :: Linux",
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    python_requires=">=3.6",
)
