
# Component versus Sample Spectral Similarity (CS3) in R

CS3 is a work-in-progress package to facilitate straightforward PARAFAC
model verification.

Code relates to the manuscript titled *Component-vs-sample similarity as
a diagnostic tool for inadequate model fitting in EEM-PARAFAC analysis.*

The current github repo is <https://github.com/MRPHarris/CS3>

### A note on code from Murphy (2011) for gradient peak detection

RamanIntegrationRange originally mentioned in this paper: Murphy, K. R.:
A Note on Determining the Extent of the Water Raman Peak in Fluorescence
Spectroscopy, Appl. Spectrosc., 65(2), 233-236, <doi:10.1366/10-06136>,
2011.

Subsequently, the function was integrated into the drEEM toolbox:
Murphy, K. R., Stedmon, C. A., Graeber, D. and Bro, R.: Fluorescence
spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 5,
6557-6566, 2013.

The function was then extracted from drEEM in MATLAB using \<open
RamanIntegrationRange.m\> Code was then ported to R, using equivalent
commands. e.g. forecast::ma() instead of smooth(). Where it was
possible, functions that preserved the ‘style’ of the code were used.
E.g. discrete numerical gradient with `pracma::gradient()`
