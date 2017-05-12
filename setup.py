
from numpy.distutils.core import Extension, setup


classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT",
    "Programming Language :: F",
    "Programming Language :: Python :: 3.5",
    "Topic :: Scientific/Engineering :: Bio-Informatics"
]


futils = Extension(
        name    = 'phd.futils.futils',
        sources = ['phd/futils/futils.f90'],
)


setup(
    name         = "phd",
    version      = "0.0.1",
    author       = "Justus Schwabedal",
    author_email = "JSchwabedal@gmail.com",
    description  = ("Python module to extract and analyze phase dynamics from oscillatory signals."),
    license      = "MIT",
    keywords     = "phase dynamics",
    url          = "https://github.com/jusjusjus/phd",
    packages     = ['phd'],
    ext_modules	 = [futils],
    classifiers	 = classifiers
)
