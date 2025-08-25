##################################################
The MultiProFit astronomical source modelling code
##################################################

DOI: `10.71929/rubin/2584108 <https://doi.org/10.71929/rubin/2584108>`_

.. abstract::

   This note describes the MultiProFit astronomical source modelling code, and particularly its development and use for modelling coadded objects in the LSST Science Pipelines.

Introduction
============

``MultiProFit`` is a software package for astronomical source modelling made by and for LSST Data Management.
This document outlines the motivation for its development, the structure of the package at the time of writing, relevant and/or unique aspects of its design, and lists some avenues for future development.
Where possible, references are made to JIRA tickets relevant to the design of the feature in question.
Comments and questions should be directed to the primary author of the package and this note, Dan S. Taranu (hereafter DST).

Motivation
==========

Many of the science cases for the Rubin Observatory/LSST survey, particularly those covered by the Galaxies and Dark Energy Science Collaborations, are dependent on measuring magnitudes (i.e. total fluxes) and shapes of galaxies.
As of 2025, the main package providing galaxy flux measurement in the science pipelines :cite:`PSTN-019` is/was  `meas_modelfit <https://github.com/lsst/meas_modelfit>`_, which dates to (at least) summer 2009.
``meas_modelfit``, like ``MultiProFit``, implements forward modelling of images of astronomical sources (star and galaxies).
Both sources and the point spread function in each image are modelled as Gaussian mixtures.
Initially, it was envisioned that the package would perform simultaneous modelling of individual exposures, dubbed "multifit" (indeed, there was a ``meas_multifit`` subpackage for some time).
By 2015, most development of ``meas_modelfit`` was completed and it had become increasingly clear that multifit would be prohibitively computationally expensive on real LSST data and provide marginal benefits over fitting coadded images.
In 2018, it was decided to focus efforts on either enhancing or replacing ``meas_modelfit`` with a package that supported simultaneous fitting of multiple bands (with one coadd per band), and with improved galaxy models.

Development of MultiProFit began in mid-2018 with the idea of extending the then relatively new `ProFit <https://github.com/ICRAR/ProFit>`_ code :cite:`2017MNRAS.466.1513R` and its fairly rudimentary Python bindings `pyprofit <https://github.com/ICRAR/pyprofit>`_.
ProFit does and did support a wider variety of parametric galaxy and PSF models, as well a variety of Bayesian optimizers, but the latter are only available in the R interface.
More importantly, while ProFit supports accurate pixel integration of various analytic profiles, it does not natively support direct evaluation of Gaussian mixture models (GMMs) - that is, all models evaluate the pixel-integrated PSF and galaxy profile, and then convolve the two via fast Fourier transforms.
Early testing confirmed that this is, broadly speaking, too slow to be feasible to run on a survey as large and deep as LSST, and GMMs are the only viable solution.
This is especially true for evaluating Sersic profiles with large indices (n>6), which necessitates oversampling of the pre-convolution profile and the PSF to compute the convolved profile accurately.

The use of GMMs in galaxy modelling dates began several years earlier with :cite:`2013PASP..125..719H`, who provided Gaussian mixture approximations to the widely-used Sersic profile; these approximations are indeed used in ``meas_modelfit``.
In 2018, the `Tractor <https://github.com/dstndstn/tractor>`_ code :cite:`2016ascl.soft04008L` was already public, and has since been used to process images from the Dark Energy Camera Legacy Survey (`DECaLS <https://www.legacysurvey.org/decamls/>`_).
However, in `DM-14637 <https://rubinobs.atlassian.net/issues/DM-14637>`_ it was determined that adopting Tractor would take at least as much effort as writing an entirely new package, and so development began on the first iteration of ``MultiProFit`` (now archived at `legacy-multiprofit <https://github.com/lsst-dm/legacy-multiprofit>`_).

Erin Sheldon's `ngmix <https://github.com/esheldon/ngmix>`_ code was also evaluated and experimentally/unofficially supported for a time through `meas_extensions_ngmix <https://github.com/lsst-dm/meas_extensions_ngmix>`_.
``ngmix``'s galaxy modelling performance was last evaluated in August 2021 in `DM-30734 <https://rubinobs.atlassian.net/issues/DM-30734>`_ and was comparable to ``multiprofit`` (at the time), with both well behind ``meas_modelfit``.
Since then, ``ngmix``'s development has focused more on uses for weak lensing such as metadetection, for which it is now included in the science pipelines; however, ``meas_extensions_ngmix`` was abandoned in favour of improving ``multiprofit``'s performance.

To our knowledge, the only other comparable GMM-based code is `forcepho <https://github.com/bd-j/forcepho>`_ :cite:`2024ascl.soft10006B`, the development of which began (unbeknownst to DST) on or prior to early 2017.
``forcepho`` has not been evaluated in any similar way, and does share some features and similarities with ``multiprofit``, as well as additional features such as support for GPU and parallel evaluation.
At this time, evaluation of its performance is neither planned nor ruled out.

Naming
======

Officially, MultiProFit is an abbreviation for Multiple Profile Fitting, calling back to its origins in ``ProFit``/``pyprofit`` (despite not having shared any code with either for many years) and indicating first-class support for multi-band, multi-exposure and multi-component fitting.
``multiprofit`` is an equally-accepted stylization and users are free to pronounce ProFit as "profit" if they prefer.
For the remainder of this note, ``multiprofit`` will be used to refer specifically to the Python package of that name, and ``MultiProFit`` more broadly to the collection of packages implementing its functionality in the Science Pipelines.

Package Structure
=================

``MultiProFit`` began as a single Python package, eventually adding C++ code with pybind11 bindings.
`DM-20193 <https://rubinobs.atlassian.net/issues/DM-20193>`_ split off the Gaussian mixture model evaluation into a separate package, `gauss2d <https://github.com/lsst-dm/gauss2d>`_, which is now a core dependency of ``multiprofit``.
`DM-30040 <https://rubinobs.atlassian.net/browse/DM-30040>`_ added `modelfit_parameters <https://github.com/lsst-dm/modelfit_parameters>`_ to provide a header-only C++ library for defining parameters with limits and transformations.
With the addition of ``modelfit_parameters``, ``gauss2d`` was further split into a extension package `gauss2d_fit <https://github.com/lsst-dm/gauss2d_fit>`_, which depends on ``modelfit_parameters`` and ``gauss2d``.
The overarching philosophy is that these packages should be standalone-installable and remain as broadly useable as possible to encourage outside contributions, though those have (unsurprisingly) not yet materialized.
Furthermore, the package division allows for ``gauss2d`` and ``gauss2d_fit`` to remain largely field-agnostic, and so astronomy-specific nomenclature and usage is left to ``multiprofit``.

To that end, ``multiprofit`` is now a standalone Python-only package, with all C++/pybind11 code in its dependencies.
``multiprofit`` does make use of `pex_config <https://github.com/lsst/pex_config>`_ and `utils <https://github.com/lsst/utils>`_ from the Science Pipelines, but these are themselves standalone packages.
Unfortunately, since ``multiprofit``'s dependencies include primarily C++ packages, it its not yet available through PyPI, but the standalone installation process is well-documented and at least sporadically tested.

``MultiProFit``'s C++ packages are provided in the Science Pipelines as third-party packages, in part because they use Meson as a build system rather than SCons, and use `doctest <https://github.com/doctest/doctest>`_ for unit tests instead of Boost (indeed, they do not use Boost at all).
These choices were in part experiments to determine the suitability of Meson and to a lesser extent doctest in new C++ packages.
DST's impressions of both are favorable, but since few new C++ packages have been written with Rubin/LSST moving closer to operations and migrating existing packages would be a considerable effort, adoption of Meson has been limited (and of doctest non-existent, because it is not a supported package or included in `rubin-env <https://anaconda.org/conda-forge/rubin-env>`_).

The last package currently in the ``MultiProFit`` family is `meas_extensions_multiprofit <https://github.com/lsst/meas_extensions_multiprofit>`_, which implements pipeline tasks following interfaces defined in `pipe_tasks <https://github.com/lsst/pipe_tasks>`_.
This name is slightly misleading, as ``MultiProFit`` is necessarily not part of the measurement plugin framework the way ``modelfit_parameters`` is - plugins are single-band only by (current) design.
However, the naming follows the convention of `meas_extensions_scarlet <https://github.com/lsst/meas_extensions_scarlet>`_, which provides interfaces for the multiband `Scarlet deblender <https://github.com/lsst/scarlet_lite>`_, itself also not a measurement plugin.

Features
========

While ``MultiProFit`` feature set has significant overlap with ``meas_modelfit`` (not to mention ``forcepho``, ``ngmix`` and ``Tractor``), this section will nonetheless list the most important features and highlight differences with alternative codes.

Gradient Evaluation
--------------------
First and foremost, the evaluation of models and their first derivatives (gradients) can be done entirely analytically - as in ``forcepho``, but in contrast with ``meas_modelfit``, where it is computed through finite differencing.
This functionality is now primarily in the ``gauss2d`` package.
Depending on the parameter values, compiler settings, etc., computing the gradients for all of the parameters of a single Gaussian component takes 2-2.5x longer than just evaluating the model.

A single Gaussian component can have up to 6 free parameters - two centroid parameters, an integral/normalization, and three shape parameters.
Thus, at most a 3x performance gain is realized when every parameter is free, and somewhat less if some of the parameters are fixed.
This is in comparison to a baseline of finite differencing, which requires one model evaluation at the new parameter values, plus at least one more for each parameter (so 7 evaluations if all parameters are free, or 5 if the centroids are fixed).

Note that ``MultiProFit`` evaluates all gradients whether or not the parameters are free, on the assumption that most will in fact be free.
This means that there is no performance benefit to fixing the centroids, though one could be realized with some additional effort.
Similarly, there's no benefit to fixing any of the ellipse parameters.
The ellipse parameter gradients share cross terms, so there is not much room for optimizing performance when only one or two of them are fixed.
In principle, some benefit could be realized if all of the ellipse parameters are fixed, but this use case is not expected to be common.

For much of ``MultiProFit``'s history, the fitting interface was implemented in Python, and so performance was actually poorer than in ``meas_modelfit``, which is implemented almost entirely in C++.
It was only after a significant effort in converting these Python classes to C++ in `DM-33219 <https://rubinobs.atlassian.net/browse/DM-33219>`_ that the overhead from Python calls was reduced enough for ``MultiProFit`` to become more performant.
Any further performance gains would require running a C++ optimizer, such as the one used in ``meas_modelfit``.
Some of the ground work for this has already been done in `DM-38617 <https://rubinobs.atlassian.net/browse/DM-38617>`_, which aims to use `GNU Scientific Library (GSL) <https://www.gnu.org/software/gsl/>`_ optimizers, but this will take considerably more effort and is not likely to be completed in 2015.

Ellipse Parameterization
------------------------
Another significant difference between ``MultiProFit`` and alternatives is in the ellipse parameterization used to represent shapes.
The widely-used `GALFIT <https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html>`_ package, for example, fits the major-axis size, the axis ratio, and a position angle, as does the original ``ProFit``.
Other fitting codes like ``meas_modelfit`` and ``ngmix`` fit the three unique elements of covariance matrix of a Gaussian, or permutations thereof more commonly used in weak lensing rather than galaxy evolution.
``MultiProFit`` is unique (to DST's knowledge) in fitting the three unique elements of the correlation matrix, i.e. an x- and y-axis size (Gaussian dispersion), and the correlation parameter ρ.
This was originally implemented in `DM-17253 <https://rubinobs.atlassian.net/browse/DM-17253>`_.
The advantage to this approach is that the correlation parameter is bounded to the range (-1, 1), and unlike the position angle, does not have periodic boundary conditions.
However, this does not remedy the general problem that shape parameter values tend to be correlated.

Another benefit of ``MultiProFit`` is that ``modelfit_parameters`` allows for the specification of a fairly wide range of parameter transformations.
The log (or log10) transformation allows limiting flux and size parameters to positive values without having to set hard limits, which many otherwise performant optimizers do not support.
Similarly, the correlation parameter ρ can be remapped to an infinite range by the logit transform.

PSF Models
----------
To fit GMMs to observed images, one must first determine a GMM representation of the PSF.
``MultiProFit`` can, of course, fit GMMs directly to images of isolated stars, but PSF models on coadded images tend to be evaluated as pixel-convolved and normalized images (representing the probability distribution function of a photon landing within a given pixel).
``MultiProFit`` allows fitting multi-component models with a single flux/normalization parameter, and a ``FluxFraction`` parameter specifying the fraction of the total remaining flux in each component.
By setting the total flux to unity, one can fit the GMM for an arbitrary combination of Gaussians, whereas in ``meas_modelfit``, the sum of component fluxes in the best-fit model is not constrained to total unity (and therefore must be re-normalized to a potentially suboptimal solution).
In practice, this method of constraining the total flux struggles with more than three components.
If any one of the flux fraction parameters tends to zero, then the remaining components in the chain also end up marginalized, and the parameter values are unavoidably correlated.
Furthermore, at present, ``MultiProFit`` is limited to supporting only two-component models in this mode, an issue to be resolved in `DM-40674 <https://rubinobs.atlassian.net/browse/DM-40674>`_

On the other hand, ``meas_modelfit`` supports higher-order Shapelet profiles for the PSF, which can model some asymmetries in the PSF.
These have yet to be implemented in ``MultiProFit``.
In practice, the values of the higher-order parameters tend to be fairly small, at least in ground-based images and when using only 2-3 components.
In the Science Pipelines, ``meas_modelfit`` is configured to fit a double Shapelet PSF, with both components essentially sharing the same axis ratio and position angle and different only in size.
``MultiProFit`` is generally configured to fit two Gaussians with independent shapes (but a shared centroid).
`DM-43357 <https://rubinobs.atlassian.net/browse/DM-43357>`_ examined the performance of both of these approaches on simulated data and found little practical difference between the two, despite varying methodologies.

MultiProFit does not yet support fitting PSF model parameters.
Fitting the PSF model simultaneously with profiles of extended objects is dangerous.
This was tested on KiDS data with `AllStarFit <https://github.com/taranu/allstarfit>`_ and the unavoidable result is that galaxies bias the PSF parameters.
It is possible to fit a common PSF model (with fractional fluxes) to multiple stars by sharing structural parameter objects.
As this methodology is not yet needed, no convenient interface for doing so is provided.
Users interested in generic PSF model fitting should either fit individual stars, or use a code like `Piff <https://github.com/lsst/meas_extensions_piff>`_ :cite:`2021ascl.soft02024J`, which, amongst other features, allows for spatial variations of PSF parameters.

MultiProFit also does not support oversampled model evaluation.
As such, it is critical to ensure that the PSF model is not undersampled - i.e., that individual Gaussian components have σ>0.8 pixels.
One way to do this when fitting stars or PSF model images is to provide a PSF model with σ>0.8.
This way, any components will be guaranteed to have sizes larger than this minimum value.

Manipulation of the PSF model parameters can allow for some flexibility in fitting PSF-convolved profiles.
When fitting a PSF-convolved model, the size parameters cannot be negative, which essentially places a hard prior on the sizes of stars.
As such, it can be helpful to slightly shrink the PSF parameters to allow stars to have smaller sizes than that of the nominal best-fit PSF.
For example, shrinking the best fit σ values by 0.01 pixels (in quadrature) means that the sizes of stars should cluster around this value, and one can determine from the bias if the PSF model size is over- or under-estimated.
The subtracted size can be added back in to the sizes of extended objects, although real galaxies are large enough that such a small offset is irrelevant, especially in quadrature.

Galaxy Profiles
---------------
More consequential than PSF modelling differences is the choice of galaxy profile.
``meas_modelfit`` implements only the fixed Sersic index profiles originally presented in :cite:`2013PASP..125..719H`.
``MultiProFit`` allows for Sersic profiles with index values between 0.5 and 6.
``meas_modelfit`` is typically configured to fit the composite model (cModel or CModel) popularized by the Sloan Digital Sky Survey (SDSS).
This model is essentially the best-fit linear combination of two independent exponential (n=1) and de Vaucouleurs (n=4) profile fits to the same galaxy.
Since linear optimization is inexpensive, this can be interpreted as a form of bulge-disk decomposition, although in practice, doing the separate nonlinear fits for each component biases the parameters of both, even for galaxies that appear composed of distinct bulge and disk components.

``MultiProFit`` allows for fitting a variety of profiles, including a single Sersic model, bulge-disk models with fixed or free Sersic indices for both parameters, and models including a central point source component.
The specific composition of the GMM representation of a given Sersic profile is not unique.
Besides varying the number of constituent Gaussians, one has flexibility in the metric used to optimize the parameters (weights and scale radii) of those parameters.

``MultiProFit``'s Sersic profile weights were originally derived in `DM-15909 <https://rubinobs.atlassian.net/issues/DM-15909>`_ and verified/validated in `DM-21287 <https://rubinobs.atlassian.net/issues/DM-21287>`_, independently from those used in ``forcepho`` and ``Tractor``.
Broadly, the weights were chosen to reproduce the 1D radial profile for a given Sersic index, excluding the very innermost and outermost regions for large values of the Sersic index.
Sersic profiles with n>6 have highly peaked central profiles, which generally do not occur in real galaxies, and which can make the profile indistinguishable from an unresolved point source.
Similarly, such profiles also have very extended outer wings, which require large number of Gaussian components to accurately model, and which can be degenerate with residual backgrounds in astronomical images (whether undersubtracted sky or any other low surface brightness feature).

MultiProFit provides two varieties of weights, one with four Gaussians and another with eight.
The eight-Gaussian variety is only necessary to reproduce the outer profile for large Sersic indices (n > 4).
Since most galaxies are faint and disk-like, by default, only the four-component version is used.

In practice, the single Sersic profile has been found to be the most robust and accurate model, even on simulated data where the galaxies are actually bulge-disk models (see e.g. `DM-42270 <https://rubinobs.atlassian.net/issues/DM-42270>`_)
Mostly this is because small and low signal-to-noise galaxies provide limited constraints on the size and shape of the smaller component (usually the bulge).
Additionally, most of the galaxies in the universe are disk-dominated and/or asymmetric, and so the addition of a bulge component adds little to the quality of the fit.
For large and/or high signal-to-noise galaxies, most optimizers tend to be sensitive to initial conditions and struggle with robustness of fits, even when additional components are statistically justifiable.

Multiband Models
----------------
``MultiProFit`` implements multiband fitting in the simplest possible fashion, by having components share structural (shape and Sersic index) parameters but fitting a separate flux normalization in each band.
The advantage to this approach is that the flux parameters remain linear and can be fit with (non-negative) least squares optimizers, which are much faster than nonlinear optimizers.

Optimizers
----------
Unlike ``forcepho``, ``ProFit`` and ``Tractor``, ``MultiProFit`` does not yet support Bayesian optimization through algorithms like MCMC.
Currently, ``MultiProFit`` uses SciPy's nonlinear maximum likelihood optimizers.
Limited and largely experimental support is also provided for `PyGMO <https://esa.github.io/pygmo2/>`_ :cite:`10.21105/joss.02338`, although none of the tested optimizers appear to perform better than SciPy's.

``meas_modelfit`` uses its own C++ optimizers, which are described in the `Doxygen C++ documentation <http://doxygen.lsst.codes/stack/doxygen/x_mainDoxyDoc/classlsst_1_1meas_1_1modelfit_1_1_optimizer.html>`_, and were inspired by Numerical Optimization :cite:`nocedal2006numerical`.

Usage
=====
The ``multiprofit`` package provides several levels of interfaces for fitting data.
These are detailed in `MultiProFit's documentation <https://pipelines.lsst.io/v/weekly/modules/lsst.multiprofit/index.html>`_.
In brief, the higher-level classes for catalog fitting implement batch fitting of individual objects given a set of images and a corresponding catalog.

``meas_extensions_multiprofit`` uses ``multiprofit``'s interfaces for batch fitting and implements subtasks for the pipeline tasks defined in ``pipe_tasks``.
The package defines methods for loading data for a single deblended object and initializing model parameters.
The tasks generate per-patch tables, with per-band tables for PSF fit parameters, and a single table for the multiband object model fits.
A subset of the columns generated by these tasks are then merged into per-patch object tables, which are then consolidated into a per-tract table (this was added in `DM-48591 <https://rubinobs.atlassian.net/browse/DM-48591>`_).

The task structure is one possible way of running multiband algorithms, defined by pipeline connections rather than through configuration settings like the single-band measurement plugin framework.
``Scarlet`` similarly implements its own multiband tasks for deblending.
Given this experience, a general multiband measurement plugin framework may not be necessary if individual stages of processing can follow a similar design of pipeline tasks as interfaces with configurable subtasks as plugins.

Future Development
==================
There are ongoing efforts to expand on and improve ``MultiProFit``'s functionality in many of the aforementioned areas.
These can be found searching for the ``multiprofit`` component `on Jira <https://rubinobs.atlassian.net/issues/?jql=project%20%3D%20%22DM%22%20AND%20component%20%3D%20%22multiprofit%22>`_, although some tickets pertain only to the dependencies (gauss2d/gauss2d_fit).

Sersic Profile weights
----------------------
The weights for the Sersic profile could be improved.
Currently, weights are defined for a fixed set of knot values and interpolated with GSL splines.
However, the knots do not yield smooth variations for every component's size and integral at all Sersic index values, which can mislead or trap optimizers using gradients.
This is to be reviewed and improved in `DM-42106 <https://rubinobs.atlassian.net/browse/DM-42106>`_

Deblending
----------
Perhaps the most anticipated upcoming addition is the re-implementation of simultaneous multi-object fitting (i.e. deblending) in `DM-42968 <https://rubinobs.atlassian.net/browse/DM-42968>`_.
Deblending - either linear, fitting only fluxes, or fully non-linear for all parameters - was implemented and tested prior to most of the performance improvements that yielded faster runtimes than ``meas_modelfit``.
Multi-object fitting in large blends is challenging - even linear deblending is memory-intensive, whereas the runtime for non-linear deblending can scale non-linearly with the number of objects.
However, this has numerous benefits, including that the outputs will not have to rely on a separate deblender like Scarlet, as long as an alternative method of parameter initialization is implemented.
Fitting models to deblended images is convenient, but there is no principled or correct way to preserve the noise in regions where neighboring objects overlap, and so it is expected that multi-object fitting will also improve parameter uncertainty estimates.

Multi-resolution fitting
------------------------
Preliminary experiments with multi-resolution (i.e. multi-survey) fitting are being conducted on `DM-46497 <https://rubinobs.atlassian.net/browse/DM-46497>`_, using HSC and HST data in COSMOS.
COSMOS data has been public for well over a decade, so this kind of joint fitting could have been implemented (in ``MultiProFit`` or other codes) many years ago.
The challenge is more practical than conceptual - it does add complexity to any code, but also requires careful curation of the input datasets.
For example, the WCS solutions for coadded images may not be consistent enough to define a common centroid in sky coordinates.
Archival data may have been processed before commonly-used star reference catalogs like Gaia were available.
Additionally, fitting stars from different epochs may require implementation of proper motions for their centroids.

References
==========

.. bibliography::