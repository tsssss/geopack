# python geopack  setup.py

from distutils.core import setup

setup(
    name= 'geopack',
    packages= ['geopack'],
    version= '1.0.1',
    description= 'Python version of geopack and Tsyganenko models',
    author= 'Sheng Tian',
    author_email= 'tianx138@umn.edu',
    url= 'https://github.com/tsssss/geopack',
    download_url= 'https://github.com/tsssss/geopack/blob/master/dist/geopack-1.0.1.tar.gz',
    requires= ['numpy','scipy'],
    platforms= ['Mac OS'],
    license= 'MIT',
    keywords= ['geopack','space physics','Tsyganenko model'],
    classifiers= [
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Operating System :: MacOS',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Physics'
    ],
    long_description= """\
Python geopack and Tsyganenko models
------------------------------------

Calculate magnetic field vector at given position and time,
based on internal models (dipole and IGRF) and external models
(Tsyganenko models, one of T89, T96, T01, T04).

Transform vectors from one coordinate to another among
GSM, GSE, GEO, GEI, SM, MAG.

Trace the given position parallel or anti-parallel to the
model magnetic field until reaching a given boundary.

This version requires Python 3 or later.
"""
)
