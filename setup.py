import glob
import os

from setuptools import setup, find_packages


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


version = 'x.y.z'
if os.path.exists('VERSION'):
    version = open('VERSION').read().strip()


setup(
    name='TransposonClassifier',
    version=version,
    description='Analyse data from RepeatModeller',
    long_description=read('README.md'),
    packages=find_packages(),
    author='Charles Nunn',
    author_email='cn14@sanger.ac.uk',
    url='https://github.com/CharlesNunn001/TransposonClassifier',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    tests_require=['nose >= 1.3'],
    install_requires=['pandas >= 0.25.3'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience  :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
    ],
)