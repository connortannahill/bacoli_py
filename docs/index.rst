.. currentmodule:: bacoli_py

Welcome
=======

:mod:`bacoli_py` is a python package for the error controlled numerical solution of 1D Parabolic PDE's. This software was built on top of modified versions of the `BACOLI <http://cs.stmarys.ca/~muir/BACOLI-3_Webpage.htm>`_ and BACOLRI Fortran packages. For any questions or comments, please email ctannahill3@gmail.com

A formal introduction to :mod:`bacoli_py` is given in the report: :download:`report <BacoliPy_Report.pdf>`.

Contents
--------

.. toctree::
    :maxdepth: 1

    tutorial
    examples/examples


Installing and Starting with bacoli_py
----------------------------------------------

Installation requires you to have the gfortran compiler (other Fortran compilers may work as well but this has not been tested). To install :mod:`bacoli_py` and all its dependencies type "pip install bacoli_py" at the command line.

:mod:`bacoli_py` is available through `PyPi <http://pypi.python.org/pypi/bacoli_py>`__. Documentation and a tutorial are available `here <tutorial.html>`__.

Source Code
-----------
The source code for :mod:`bacoli_py` can be found on `GitHub <https://github.com/connortannahill/bacoli_py>`__.

Compiling on Windows
----------------------------------
To install on windows:

#. Get MinGW with gfortran `here <http://www.equation.com/servlet/equation.cmd?fa=fortran>`__.

#. Compile from `source <https://github.com/connortannahill/bacoli_py>`__ using ``python setup.py config --compiler=mingw32 build --compiler=mingw32 install``

Bug Reports
-----------
Please report any issues you experience while using this software on GitHub `here <https://github.com/connortannahill/bacoli_py/issues>`__.
