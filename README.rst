Decomp
======

Decompose the optical emissions from disk and jet with reveberation mapping data.

An application to 3C 273 was published in 
`Li, Y.-R., et al. 2020, ApJ, 897, 18 <https://ui.adsabs.harvard.edu/abs/2020ApJ...897...18L/abstract>`_.

Usage
-----

First prepare the light curve data for the driving continuum, emission line, and radio.
Then edit the parameter file ``src/param`` to specify the file names for these data, as well as 
to set the anticipated range of time lags between the continumm and emision line, and between the continuum 
and radio.

Run the code with the following command,

.. code-block:: bash

    mpiexec -n np ./decomp src/param

Here, ``np`` is the number of cores used for running. Replace it with a number you want to use. 