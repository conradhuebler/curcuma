#!/usr/bin/env python3
"""
Setup script for Curcuma Python bindings

This setup.py uses scikit-build-core to build the C++ extension module
via CMake and pybind11.

Installation:
    pip install .                 # Standard installation
    pip install -e .              # Editable/development installation
    pip install -v .              # Verbose build output

Requirements:
    - Python >= 3.7
    - CMake >= 3.18
    - C++14 or C++17 compiler
    - pybind11 >= 2.11

Claude Generated: Python package setup for Curcuma
Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
"""

from setuptools import setup, find_packages
import os
import sys

# Read version from git or use default
try:
    import subprocess
    git_tag = subprocess.check_output(
        ["git", "describe", "--tags"],
        cwd=os.path.dirname(__file__) or ".",
        stderr=subprocess.DEVNULL
    ).decode().strip()
    __version__ = git_tag.lstrip("v")
except:
    __version__ = "0.1.0.dev"

# Read README for long description
readme_path = os.path.join(os.path.dirname(__file__), "..", "README.md")
if os.path.exists(readme_path):
    with open(readme_path, "r", encoding="utf-8") as f:
        long_description = f.read()
else:
    long_description = "Curcuma - Molecular Modelling and Simulation Toolkit"

setup(
    name="curcuma",
    version=__version__,
    author="Conrad Hübler",
    author_email="Conrad.Huebler@gmx.net",
    description="Python interface to the Curcuma computational chemistry toolkit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/conradhuebler/curcuma",
    license="GPL-3.0",

    # Python package configuration
    packages=find_packages(where="python"),
    package_dir={"": "python"},

    # Requirements
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.19",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "flake8",
        ],
        "docs": [
            "sphinx>=4.0",
            "sphinx-rtd-theme",
        ],
    },

    # Classifiers
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: C++",
    ],

    # Entry points
    entry_points={
        "console_scripts": [
            "curcuma-python=curcuma.cli:main",
        ],
    },

    # Include data files
    include_package_data=True,
    zip_safe=False,
)
