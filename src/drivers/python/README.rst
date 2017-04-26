FABM and python
=======================

This folder contains source and build tools to compile and build FABM as a
python module. In addition a number of useful scripts are included.

----

When configuring for compilation using `CMake <http://cmake.org>`_ and additional
target 'wheel' is added - either in VisualStudio or in Makefile - depending
on platform. The taget 'wheel' will build a '.whl' file that can be installed
using `pip <https://pip.pypa.io/en/stable/>`_.

This is mainly for developers/release managers as a few additional steps are still required.

It is the intention to have up-to-date wheels for some platforms on `release <https://github.com/fabm-model/fabm/releases>`_.
