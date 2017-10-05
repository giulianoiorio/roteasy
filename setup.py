from setuptools import setup
import numpy
import shutil




setup(
    name='roteasy',
    version='1.6.0',
    author='Giuliano Iorio',
    author_email='giuliano.iorio@unibo.it',
    url='https://github.com/iogiul/roteasy.git',
    packages=['roteasy','roteasy/src'],
    install_requires=['numpy>=1.9',],
    include_dirs=[numpy.get_include()]
)


shutil.rmtree('build')
shutil.rmtree('dist')
shutil.rmtree('roteasy.egg-info')
