Development
***********

If you are interested in contributing to the PyPGx package, please follow the instructions listed below.

First, you need to clone the GitHub repository for PyPGx.

.. code-block:: console

   git clone https://github.com/sbslee/pypgx

Move to the ``pypgx`` directory and then create a new branch called ``my_branch``. Checkout the branch.

.. code-block:: console

   cd pypgx
   git branch my_branch
   git checkout my_branch

Create a fresh conda environment called ``my_env``. Make sure to activate the environment.

.. code-block:: console

   conda create -n my_env --file requirements.txt
   conda activate my_env

Next, install PyPGx in the development mode.

.. code-block:: console

   pip install -e .

Now, any changes you make in the ``pypgx`` directory will be automatically reflected. Once you are finished, run the unit test to make sure everything still works as expected.

.. code-block:: console

   pytest

If and only if the test is passed, make a pull request and it will be reviewed as soon as possible.
