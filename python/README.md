The libcasm-mapping package
===========================

The `libcasm-mapping` Python package provides a Python interface to the CASM structure mapping library.

This version of `libcasm-mapping` is compatible with version 2.X of [`CASMcode_mapping`](https://github.com/prisms-center/CASMcode_mapping/).


Install
=======

Installation of `libcasm-mapping` requires:
- Python >=3.8
- The compatible version of the CASM C++ mapping library is already installed.
- Development environment that allows compiling the pybind11 interface to CASM C++ (i.e. C++ compiler with support for c++17)

Normal installation:

    pip install .

Editable installation:

    pip install -e .


Building documentation
======================

Install documentation requirements:

    pip install -r doc_requirements.txt

Install `libcasm-mapping`

Build and open the documentation:

    cd doc
    make html
    open _build/html/index.html


Testing
=======

To install testing requirements, do:

    pip install -r test_requirements.txt

Use `pytest` to run the tests. To run all tests, do:

    python -m pytest -rsap tests

As an example of running a specific test, do:

    python -m pytest -rsap tests/test_prim.py::test_asymmetric_unit_indices


Formatting
==========

Use yapf. From CASMcode_mapping/python do:

    yapf -ir .
