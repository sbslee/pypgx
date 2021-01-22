Development
***********

If you are interested in contributing to the pypgx package, please follow the instructions listed below.

First, you need to clone the GitHub repository for pypgx.

.. code-block::

   git clone https://github.com/sbslee/pypgx

Move to the ``pypgx`` directory and then create a new branch called ``my_branch``. Checkout the branch.

.. code-block::

   cd pypgx
   git branch my_branch
   git checkout my_branch

Create a fresh conda environment called ``my_env``. Make sure to activate the environment.

.. code-block::

   conda create -n my_env --file requirements.txt
   conda activate my_env

Next, install pypgx in the development mode.

.. code-block::

   pip install -e .

Now, any changes you make in the ``pypgx`` directory will be automatically reflected. Once you are finished, run the unit test to make sure everything still works as expected.

.. code-block::

   pytest

If and only if the test is passed, make a pull request and it will be reviewed as soon as possible.
