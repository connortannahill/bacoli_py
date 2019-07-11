# Authors: Connor Tannahill <ctannahill3@gmail.com> and Paul Muir, 2019

from setuptools import setup, Extension
# from distutils.core import setup
from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import Configuration
from sphinx.setup_command import BuildDoc as SphinxBuildDoc


descr = """
bacoli_py is a Python package for solving systems of 1D Parabolic
Partial Differential Equations. Wraps a slightly modified version of BACOLI
(see http://cs.stmarys.ca/~muir/BACOLI-3_Webpage.htm) and BACOLRI. Tutorial,
documentation and examples can be found at https://bacoli-py.readthedocs.io/en/latest/.

For any questions or suggestions please email Connor Tannahill at
ctannahill3@gmail.com
"""

VERSION             = '1.0'
DISTNAME            = 'bacoli_py'
DESCRIPTION         = 'Python Package for Solving 1D PDEs'
LONG_DESCRIPTION    = descr
URL                 = 'https://pypi.org/project/bacoli_py/'
LICENSE             = 'BSD License'
AUTHOR              = 'Connor Tannahill and Paul Muir'
AUTHOR_EMAIL        = 'ctannahill3@gmail.com'
DOWNLOAD_URL        = URL
EXTRA_INFO          = dict(
    install_requires=['numpy'],
    classifiers=['Development Status :: 5 - Production/Stable',
                 'Programming Language :: Python',
                 'License :: OSI Approved :: BSD License',
                 'Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering',
                 'Topic :: Scientific/Engineering :: Mathematics',
                 'Operating System :: OS Independent']
)

def configuration(parent_package='', top_path=None):
    config = Configuration(DISTNAME, parent_package, top_path, 
                           version = VERSION,
                           author = AUTHOR,
                           author_email = AUTHOR_EMAIL,
                           description = DESCRIPTION,
                           license = LICENSE,
                           url = URL)
                           # long_description = LONG_DESCRIPTION,
                           # long_description_content_type='text/plain')

    config.add_data_files('LICENSE.txt')
    config.add_data_files('bacoli_py/lib/BACOLI-LICENSE.txt')
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
    version='1.0',
    license='needed',
    long_description=open('README.txt').read(),
    ext_modules=[bacoli_interface]
)