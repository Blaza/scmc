SCMC - Stochastic Collocation Monte Carlo
================

This is a development version of the `scmc` package for the R
programming language.

## Installation

The package isn’t available on CRAN so the only way to install the
package is to use the `devtools` package and run

``` r
devtools::install_github("blaza/scmc")
```

## Usage

The main function currently implemented is `univariate_sampler` which is
a flexible implementation of the method (and thus with a bit more
complicated interface) for generating univariate distributions. We’ll
cover here a couple of basic examples which give an overall picture of
the package capabilities.

### Example: Logistic distribution

We’ll generate variates from the [Logistic
distribution](https://en.wikipedia.org/wiki/Logistic_distribution).

The SCMC method implies interpolating <img src="/tex/4977c5a48089ee730209022980df931a.svg?invert_in_darkmode&sanitize=true" align=middle width=118.26284084999999pt height=26.76175259999998pt/>, where <img src="/tex/91aac9730317276af725abd8cef04ca9.svg?invert_in_darkmode&sanitize=true" align=middle width=13.19638649999999pt height=22.465723500000017pt/> is
the target random variable and <img src="/tex/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode&sanitize=true" align=middle width=14.908688849999992pt height=22.465723500000017pt/> is a random variable which can be
efficiently generated, and generating the samples <img src="/tex/1e84c29c315b23141c448eb0f38c269e.svg?invert_in_darkmode&sanitize=true" align=middle width=93.29428679999998pt height=21.68300969999999pt/> using
the formula <img src="/tex/e2c91384c1e794cae27f5526ed057013.svg?invert_in_darkmode&sanitize=true" align=middle width=162.59677829999998pt height=26.76175259999998pt/>, where <img src="/tex/1b343a3f3b853d91167a5af80932c8be.svg?invert_in_darkmode&sanitize=true" align=middle width=5.663225699999993pt height=21.68300969999999pt/> are variates from
the <img src="/tex/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode&sanitize=true" align=middle width=14.908688849999992pt height=22.465723500000017pt/> distribution. In this example we’ll use the standard normal
variable <img src="/tex/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode&sanitize=true" align=middle width=14.908688849999992pt height=22.465723500000017pt/>. By default, we use the `RcppZiggurat::zrnorm` function to
generate normal variates.

The code to generate the Logistic distribution in the `scmc` package is

``` r
library(scmc)
# create the sampler
sampler <- univariate_sampler(qlogis, gaussian_nodes(7))
```

    ## Loading required package: RcppZiggurat

``` r
# generate 10000 random variates
smp <- sampler(1e5)
```

In its basic form, the `univariate_sampler` function requires the
inverse <img src="/tex/926c661e425ad32f4230f7859416fdba.svg?invert_in_darkmode&sanitize=true" align=middle width=48.79468769999999pt height=26.76175259999998pt/> (i.e. the quantile function of <img src="/tex/91aac9730317276af725abd8cef04ca9.svg?invert_in_darkmode&sanitize=true" align=middle width=13.19638649999999pt height=22.465723500000017pt/>) as the first
argument, and the nodes for the interpolation. In cases where normally
distributed <img src="/tex/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode&sanitize=true" align=middle width=14.908688849999992pt height=22.465723500000017pt/> are used, optimal nodes for interpolation are the nodes
of the Gaussian quadrature with respect to the weight function <img src="/tex/510e72fbfa14c5f8620aa97e8cb52288.svg?invert_in_darkmode&sanitize=true" align=middle width=52.82431934999998pt height=24.65753399999998pt/>
(density of <img src="/tex/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode&sanitize=true" align=middle width=14.908688849999992pt height=22.465723500000017pt/>). The third argument to `univariate_sampler` is `xdist`
which is by default `"norm"`, indicating the standard normal
distribution.

The quality of the generated sample can be visualized with it’s density

``` r
# plot the sample density
plot(density(smp))

# add a plot of the theoretical logistic density
curve(dlogis(x), add = TRUE, col = "green")
```

![](README.tex_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

The curves are nearly the same, so the approximation is good.
