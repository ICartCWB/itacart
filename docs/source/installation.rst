============
Installation
============

This guide provides instructions on how to install the `itacart` package.

Using pip
---------

The recommended way to install `itacart` is from the Python Package Index (PyPI) using pip.

.. code-block:: bash

   pip install itacart

This will download and install the latest stable version of the package and its required dependencies.

From Source (for Developers)
----------------------------

If you want to contribute to the project or need the latest development version, you can install it directly from the source code on GitHub.

1. **Clone the repository:**

   .. code-block:: bash

      git clone https://github.com/ITACaRT/itacart.git
      cd itacart

2. **Install in editable mode:**

   .. code-block:: bash

      pip install -e .

   The `-e` or `--editable` flag allows you to modify the source code and have the changes immediately reflected in your installed package without needing to reinstall.

Verify Installation
-------------------

To verify that the installation was successful, you can run the following command in your Python interpreter:

.. code-block:: python

   import itacart
   print(itacart.__version__)

This should print the installed version of the package.
