
# Component versus Sample Spectral Similarity (CS3) in R

Current main branch version: 0.1.1

CS3 is a work-in-progress package to facilitate straightforward PARAFAC
model verification, building on the R fluorescence framework provided by
the EEM, eemR, and staRdom packages.

Code relates to the manuscript titled *Component-vs-sample similarity as
a diagnostic tool for inadequate model fitting in EEM-PARAFAC analysis.*

The current github repo is <https://github.com/MRPHarris/CS3>.

### Package framework

A single function is exported from this package, `per_eem_ssc`. This
function compares a PARAFAC component with constituent EEM spectra at
the target component peak wavelength position. Four corrections are
optionally applied during this operation:

1)  Spectral correction. During spectral correction, the contribution
    from non-target components is subtracted from the sample spectra
    being compared. In this manner, the target component is compared
    only against sample + residual fluorescence, and similarity
    differences are not contributed to by other components.

2)  Interpolation. Spectra are interpolated to an evenly spaced 1nm
    bandwidth.

3)  Savitzky-Golay smoothing. Raw sample spectra are noisy - part of the
    strength of PARAFAC is that this noise can be removed. Thus in
    directly comparing PARAFAC and sample spectra, it is often necessary
    to smooth the data to some degree. SG smoothing preserves the shape
    of the spectra whilst dispensing with a lot of noisiness.

4)  Incomplete excitation peak removal. The shift and shape sensitive
    congruence metric was developed by Wunsch et al., 2019, in order to
    better compare excitation spectra. However, this metric is also
    useful in picking up small excitation spectral differences. To apply
    it, it is necessary to remove incomplete peaks so as not to bias the
    peak position penalty term (incorrect assignment of peaks). A
    gradient-based peak detection and trimming method based upon
    Murphy (2011) is applied for this purpose. SSC comparisons are made
    only to trimmed spectra.

After passing through some combination of these corrections, spectra are
evaluated for TCC and/or SSC, and the results returned. High metric
values indicate that a component strongly matches the spectra in the
EEMs it claims it is explaining. Low values indicate that this component
is poorly fitted. However, the values should not be averaged and taken
at face value, as this would essentially be a sub-par split-half test.
Instead, the values should be plotted against other model performance
metrics along with the component loadings, and used to infer where the
model is falling down. Are there specific samples where the model and
sample similarity falls apart? Why might this be? Is spectral averaging
occurring; are more components required?

### Using the package

The `per_eem_ssc` only requires two things: a PARAFAC model object,
returned from `staRdom::eem_parafac()`, and the eemlist used to generate
the model. The eemlist must be compliant with the staRdom/EEM/eemR
fluorescence framework for R (see Pucher et al., 2019). Supply these
parameters to the function, and you’re good to go. It is recommended
that the function be run repeatedly with and without correction
parameters in order to compare the results.

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

#### References

Murphy, K. R. (2011). A Note on Determining the Extent of the Water
Raman Peak in Fluorescence Spectroscopy. Applied Spectroscopy, 65(2),
233–236. <https://doi.org/10.1366/10-06136>

Pucher, M., Wünsch, U., Weigelhofer, G., Murphy, K., Hein, T., &
Graeber, D. (2019). staRdom: Versatile Software for Analyzing
Spectroscopic Data of Dissolved Organic Matter in R. Water, 11, 2366.
<https://doi.org/10.3390/w11112366>

Wünsch, U. J., Bro, R., Stedmon, C. A., Wenig, P., & Murphy, K. R.
(2019). Emerging patterns in the global distribution of dissolved
organic matter fluorescence. Analytical Methods, 11(7), 888–893.
<https://doi.org/10.1039/C8AY02422G>
