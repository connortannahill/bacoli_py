# Authors: Connor Tannahill <ctannahill3@gmail.com> and Paul Muir, 2019

from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import Configuration
from sphinx.setup_command import BuildDoc as SphinxBuildDoc

descr = """\
bacoli_py is a Python package for solving systems of 1D Parabolic
Partial Differential Equations. Wraps a slightly modified version of BACOLI
(see http://cs.stmarys.ca/~muir/BACOLI-3_Webpage.htm) and BACOLRI. Tutorial,
documentation and examples can be found at http://pypi.python.org/pypi/bacoli_py.

For any questions, suggestions or suggestions please email Connor Tannahill at
ctannahill3@gmail.com
"""

VERSION             = '1.0'
DISTNAME            = 'bacoli_py'
DESCRIPTION         = 'Python Package for Solving 1D PDEs'
LONG_DESCRIPTION    = descr
MAINTAINER          = 'Connor Tannahill',
MAINTAINER_EMAIL    = 'ctannahill3@gmail.com',
URL                 = '**Needed'
LICENSE             = 'BSD Liscense'
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
                           packages = ['bacoli_py'],
                           version = VERSION,
                           maintainer = MAINTAINER,
                           maintainer_email = MAINTAINER_EMAIL,
                           description = DESCRIPTION,
                           license = LICENSE,
                           url = URL,
                           long_description = LONG_DESCRIPTION)
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
    name='bacoli_py',
    version='1.0dev',
    packages=['bacoli_py'],
    license='needed',
    long_description=open('README.txt').read(),
    ext_modules=[bacoli_interface]
)