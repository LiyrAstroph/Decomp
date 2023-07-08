Decomp
======

Decompose the optical emissions from disk and jet with reveberation mapping data.

An application to 3C 273 was published in 
`Li, Y.-R., et al. 2020, ApJ, 897, 18 <https://ui.adsabs.harvard.edu/abs/2020ApJ...897...18L/abstract>`_.

Installation
------------
Third-party packages required:  

* MPICH
* Lapack/Lapacke
* GSL
* BLAS
* CDNest (https://github.com/LiyrAstroph/CDNest)

After installing the above packages, edit the Makefile to set appropriate configurations to the package paths. 
Then install ``decomp`` using the command 

.. code-block:: bash

  make 

This will create an excutable binary file ``decomp``.


Usage
-----

First prepare the light curve data for the driving continuum, emission line, and radio.
Then edit the parameter file ``src/param`` to specify the file names for these data, as well as 
to set the anticipated range of time lags between the continumm and emision line, and between the continuum 
and radio.  Place these files into the subdirectory ``data/``.

Run the code with the following command,

.. code-block:: bash

    mpiexec -n np ./decomp src/param

Here, ``np`` is the number of cores used for running. Replace it with a number you want to use. 

The main outputs of the code include

* ``posterior_sample.txt``： posterior sample of parameters. 
  
  The column orders of parameters are

  0-2：systematic errors of continuum, line, and radio

  3-4: continuum DRW varaibility parameters 

  5-6: radio DRW varaibility parameters 

  7-9: line transfer function (Gaussian) parameters, amplitude, center, sigma

  10-12: radio tranfer function (Gaussian) parameters, amplitude, center, sigma

  *Note that center can be used to indicate the time lag.*

* ``con_rec.txt``, ``cond_rec.txt``, ``conj_rec.txt``: reconstrunction to total, disk and jet component
  of the continuum light curve.

* ``line_rec.txt``: reconstrunction to the line light curve.
  
* ``radio_rec.txt``: reconstrunction to the radio light curve.
   