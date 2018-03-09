---
title: WeibullR
---

An R package for Life Data Analysis
-----------------------------------

 

The WeibullR package provides a flexible data entry capability with three levels
of usage.

 

### *Quick Fit* Functions

Functions with intuitive names `MLEw2p` through `MRRln3p` for preparing simple
fits, bounds, and displays using default options. Only data sets with exact
failure times and/or suspensions are processed.

The quick fit functions return a simple named vector of the fitted parameters
with appropriate goodness of fit measure(s).

Optional preparation of appropriate interval bounds (at 90\\% confidence), or a
display of fit and bounds are controlled by two final arguments taking logical
entry, such that a function call like `MLEw2p(input_data,T,T)` will generate a
plot with the fitted data and confidence interval bounds.

When the first logical for bounds is set to TRUE, the returned object will be a
list with the fitted parameter vector first and dataframe of bound values
second.

 

### wblr Object Model

Construction of a wblr object is initiated by providing a data set through
function `wblr`.

Modification of the object with the progressive addition of fits and confidence
interval bounds is made via functions `wblr.fit` and `wblr.conf`.

Fine control over many aspects of fit, confidence, and display are made possible
using a flexible options mechanism.

Display for single object models is via S3 methods `plot` or `contour`, while
multiple objects\* (provided as a list) \*can be displayed on a single plot
using`plot.wblr`, `plot_contour`, or `contour.wblr`.

 

### Back-end Functions

Access to back-end functions providing all the functionality of the upper levels
of usage are provided as exported functions.

These functions may provide advanced users with resources to expand analysis
further than has been implemented within the WeibullR package.

 

Data Entry
----------

Data entry is made through the *Quick Fit* functions, `wblr`, or on the backend
through `getPPP` for rank regression, and `mleframe` for mle processing.

In all cases the primary argument `x` can be a vector of exact time failures or
a dataframe with `time`, and `event`columns as a minimum. An additional column
`qty` may optionally be used to record duplicated data.

If the dataframe entry is not used (in favor of an exact time failure vector), a
second argument `s`can be used to enter a vector of last observed success times
for right censored data (suspensions).

Beyond the entry of the first two data types, interval data (including
discoveries with last known success time=0) are entered via argument `interval`
as a dataframe with columns`left`, and `right` as a miniumum.

As with the primary argument dataframe entry, an additional column `qty` may
optionally be used to record duplicated interval data.

Such interval data entry is not supported with the Quick Fit functions.
