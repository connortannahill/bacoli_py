# Authors: Connor Tannahill <ctannahill3@gmail.com> and Paul Muir, 2019

from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import Configuration

descr = """
bacoli_py is a Python package for the error controlled numerical solution of 1D time-dependent PDEs. wraps a slightly modified version of the BACOLI
(see http://cs.stmarys.ca/~muir/bacoli-3_webpage.htm) and BACOLRI solvers, written in Fortran. a tutorial,
documentation, and examples can be found at https://bacoli-py.readthedocs.io/en/latest/.

for any questions or suggestions please email Connor Tannahill at
ctannahill3@gmail.com
"""

version             = '1.0'
distname            = 'bacoli_py'
description         = 'Python package for the error controlled numerical solution to 1D time-dependent PDEs'
long_description    = descr
url                 = 'https://pypi.org/project/bacoli_py/'
license             = 'BSD license'
author              = 'Connor Tannahill and Paul Muir'
author_email        = 'ctannahill3@gmail.com'
download_url        = url
classifiers=['Development Status :: 5 - Production/Stable',
             'Programming Language :: Python',
             'Intended Audience :: Science/Research',
             'Topic :: Scientific/Engineering',
             'Topic :: Scientific/Engineering :: Mathematics',
             'Operating System :: OS Independent']

def configuration(parent_package='', top_path=None):
    config = Configuration(distname, parent_package, top_path, 
                           packages=['bacoli_py'],
                           version = version,
                           author = author,
                           author_email = author_email,
                           description = description,
                           license = license,
                           url = url)
                           # long_description = long_description,
                           # long_description_content_type='text/plain')


    config.add_data_files('LICENSE.txt')
    config.add_data_files('bacoli_py/lib/BACOLRI-LICENSE.txt')
    config.add_data_files('bacoli_py/lib/BACOLRI-LICENSE.txt')

    return config

bacoli_interface = Extension('bacoli_interface',
                         sources=['bacoli_py/lib/bacoli_interface.f95',
                                  'bacoli_py/lib/bacoli.f',
                                  'bacoli_py/lib/bacoli-aux.f',
                                  'bacoli_py/lib/d1mach_i1mach.f95',
                                  'bacoli_py/lib/bacolri.f',
                                  'bacoli_py/lib/bacolri-aux.f',
                                  'bacoli_py/lib/bacoli_interface.pyf'])


setup(configuration=configuration,
    version=version,
    # packages=setuptools.find_packages(),
    include_package_data = True,
    long_description=open('README.txt').read(),
    requires=['numpy'],
    classifiers=classifiers,
    platforms=['any'],
    ext_modules=[bacoli_interface]
)