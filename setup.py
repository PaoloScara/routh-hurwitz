from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="routh-hurwitz",
    version="2.0.0",
    author="Paolo Scaramuzza",
    description="Routh-Hurwitz stability criterion analyzer with symbolic computation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/PaoloScara/routh-hurwitz",
    py_modules=["routh"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.8",
    install_requires=[
        "sympy>=1.12",
        "numpy>=1.24.0",
    ],
)
