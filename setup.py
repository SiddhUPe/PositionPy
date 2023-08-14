"""A setuptools based setup module.

See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
   
    name="PositionPy",
    version="1.0.0",  # Required
    description="PositionPy is a simple python library geared for User Friendly calculations for the positions of the Sun, Planets, Stars, and more.",  # Optional
  
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Siddhu Pendyala",
    author_email="v.siddhu.pendyala@gmail.com",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Build Tools",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3 :: Only",
    ],
    
    keywords="Astronomy, Sun, Science, Calculation, Ephemeris",
    packages=find_packages(where="PositionPy"),
    python_requires=">=3.7, <4",
    install_requires=["astroquery", "Solarflare"],
)
