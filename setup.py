# python geopack  setup.py
import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='geopack',
    version='1.0.7',
    author='Sheng Tian',
    author_email='tianx138@umn.edu',
    description='Python version of geopack and Tsyganenko models, compatible with geopack05 and geopack08',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url= 'https://github.com/tsssss/geopack',
    requires= ['numpy','scipy'],
    platforms= ['Mac OS'],
    license= 'MIT',
    keywords= ['geopack','space physics','Tsyganenko model'],
    packages= setuptools.find_packages(),
    package_data={'':['*.txt','*.md']},
    classifiers= [
        'Programming Language :: Python :: 3',
        'Operating System :: MacOS',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Physics'
    ],
)
