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

The SCMC method implies interpolating <img src="/tex/e2e4d72726b98b6c000ee11c5a2083be.svg?invert_in_darkmode&sanitize=true" align=middle width=88.53565544999998pt height=28.894955100000008pt/>, where <img src="/tex/91aac9730317276af725abd8cef04ca9.svg?invert_in_darkmode&sanitize=true" align=middle width=13.19638649999999pt height=22.465723500000017pt/>
is the target random variable and <img src="/tex/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode&sanitize=true" align=middle width=14.908688849999992pt height=22.465723500000017pt/> is a random variable which can
be efficiently generated, and generating the samples
<img src="/tex/51688ec4a2690fd0ad97e4e22eca84aa.svg?invert_in_darkmode&sanitize=true" align=middle width=108.51391319999999pt height=21.68300969999999pt/> using the formula <img src="/tex/bd85552af787a7f853c8d41fa91ca8ae.svg?invert_in_darkmode&sanitize=true" align=middle width=127.25513459999999pt height=28.894955100000008pt/>,
where <img src="/tex/5add1d368d6bcc924b8b5b96abe9b68e.svg?invert_in_darkmode&sanitize=true" align=middle width=11.84271164999999pt height=22.831056599999986pt/> are variates from the <img src="/tex/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode&sanitize=true" align=middle width=14.908688849999992pt height=22.465723500000017pt/> distribution. In this
example we’ll use the standard normal variable <img src="/tex/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode&sanitize=true" align=middle width=14.908688849999992pt height=22.465723500000017pt/>. By default, we use
the `RcppZiggurat::zrnorm` function to generate normal variates.

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
inverse <img src="/tex/68f61902b2e709d867f41bf7747c4e7c.svg?invert_in_darkmode&sanitize=true" align=middle width=29.680490399999993pt height=28.894955100000008pt/> (i.e. the quantile function of <img src="/tex/91aac9730317276af725abd8cef04ca9.svg?invert_in_darkmode&sanitize=true" align=middle width=13.19638649999999pt height=22.465723500000017pt/>) as the first
argument, and the nodes for the interpolation. In cases where normally
distributed <img src="/tex/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode&sanitize=true" align=middle width=14.908688849999992pt height=22.465723500000017pt/> are used, optimal nodes for interpolation are the
nodes of the Gaussian quadrature with respect to the weight function
<img src="/tex/c1de371a4982be00c4ba31c772465407.svg?invert_in_darkmode&sanitize=true" align=middle width=42.72503069999999pt height=24.65753399999998pt/> (density of <img src="/tex/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode&sanitize=true" align=middle width=14.908688849999992pt height=22.465723500000017pt/>). The third argument to
`univariate_sampler` is `xdist` which is by default `"norm"`, indicating
the standard normal distribution.

The quality of the generated sample can be visualized with it’s density

``` r
# plot the sample density
plot(density(smp))

# add a plot of the theoretical logistic density
curve(dlogis(x), add = TRUE, col = "green")
```

![](README.tex_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

The curves are nearly the same, so the approximation is good.

### Example: Gamma distribution

For the next example, we’ll use the [gamma
distribution](https://en.wikipedia.org/wiki/Gamma_distribution),
specifically <img src="/tex/dc497f9bd34f58ba60395b2ad6875709.svg?invert_in_darkmode&sanitize=true" align=middle width=46.80373994999999pt height=24.65753399999998pt/>. This is a positive distribution, so we
would like to transform it using the `log` transform to get a real
valued random variable and then upon sampling use the `exp` transform to
get a sample from the original distribution. Basically, we model
<img src="/tex/816ffca0f61d00f9b149e524bf3e07a6.svg?invert_in_darkmode&sanitize=true" align=middle width=37.16895599999999pt height=22.831056599999986pt/> instead of <img src="/tex/91aac9730317276af725abd8cef04ca9.svg?invert_in_darkmode&sanitize=true" align=middle width=13.19638649999999pt height=22.465723500000017pt/> and sample <img src="/tex/310638b5e7d0326481035b4da99224b2.svg?invert_in_darkmode&sanitize=true" align=middle width=37.46096969999999pt height=27.91243950000002pt/> to get <img src="/tex/91aac9730317276af725abd8cef04ca9.svg?invert_in_darkmode&sanitize=true" align=middle width=13.19638649999999pt height=22.465723500000017pt/>. The
code example follows

``` r
library(scmc)
# create the sampler
sampler <- univariate_sampler(function(x) qgamma(x, 5, 2),
                              gaussian_nodes(7),
                              transform = log, # the transformation of
                                               # the quantile function
                              inv_transform = exp) # the inverse of transform

# generate 10000 random variates
smp <- sampler(1e5)
```

We can check the quality of the sample distribution by plotting the
empirical cdf and the theoretical cdf of <img src="/tex/91aac9730317276af725abd8cef04ca9.svg?invert_in_darkmode&sanitize=true" align=middle width=13.19638649999999pt height=22.465723500000017pt/> or the histogram
overlayed with the density.

``` r
par(mfrow=c(1, 2))
smp_ecdf <- ecdf(smp)
curve(smp_ecdf(x), xlim = c(0, 7))
curve(pgamma(x, 5, 2), add = TRUE, col = "green")

hist(smp, breaks = 100, probability = TRUE)
curve(dgamma(x, 5, 2), add = TRUE, col = "green")
```

![](README.tex_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Again, the approximation is excellent.

### Example: Students t distribution

Now we’ll deomnstrate using the grid stretching technique for
heavy-tailed distributions. We use the [students t
distribution](https://en.wikipedia.org/wiki/Student%27s_t-distribution)
with 2 degrees of freedom. All that’s needed is to add a `gss` argument
which specifies the <img src="/tex/8cda31ed38c6d59d14ebefa440099572.svg?invert_in_darkmode&sanitize=true" align=middle width=9.98290094999999pt height=14.15524440000002pt/> value in the technique (section 4.1.1. in
Grzelak et al. 2014).

``` r
library(scmc)
# create the sampler
sampler <- univariate_sampler(function(x) qt(x, df = 2),
                              gaussian_nodes(15),
                              gss = 1.657)

# generate 10000 random variates
smp <- sampler(1e5)
```

We’ll use the Kolmogorov-Smirnov test now to test the distribution

``` r
ks.test(smp, function(x) pt(x, df = 2))
```

    ## 
    ##  One-sample Kolmogorov-Smirnov test
    ## 
    ## data:  smp
    ## D = 0.0024044, p-value = 0.6098
    ## alternative hypothesis: two-sided

A large p-value indicates a good distribution.

### Example: Beta distribution

This time we demonstrate sampling from a bounded distribution. A good
example is the [beta
distribution](https://en.wikipedia.org/wiki/Beta_distribution). It comes
in a wide variety of shapes, and we’ll use the <img src="/tex/2b2005e0d5ff12d4bf181f11169529b2.svg?invert_in_darkmode&sanitize=true" align=middle width=75.39400604999999pt height=24.65753399999998pt/>
distribution, with a distinct “U” shape. For this example we use the
chebyshev points for interpolation, and interpolate directly
<img src="/tex/68f61902b2e709d867f41bf7747c4e7c.svg?invert_in_darkmode&sanitize=true" align=middle width=29.680490399999993pt height=28.894955100000008pt/>. This is equivalent to using a uniformly distributed <img src="/tex/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode&sanitize=true" align=middle width=14.908688849999992pt height=22.465723500000017pt/>,
so we’ll set the `xdist="unif"` argument.

``` r
library(scmc)
# create the sampler
sampler <- univariate_sampler(function(x) qbeta(x,  0.5, 0.5),
                              chebyshev_nodes(15),
                              xdist = "unif")

# generate 10000 random variates
smp <- sampler(1e5)
```

And we visualise the distribution with the histogram and run `ks.test`.

``` r
hist(smp, breaks = 100, probability = TRUE)
curve(dbeta(x, 0.5, 0.5), add = TRUE, col = "green")
```

![](README.tex_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
ks.test(smp, function(x) pbeta(x, 0.5, 0.5))
```

    ## Warning in ks.test(smp, function(x) pbeta(x, 0.5, 0.5)): ties should not be
    ## present for the Kolmogorov-Smirnov test

    ## 
    ##  One-sample Kolmogorov-Smirnov test
    ## 
    ## data:  smp
    ## D = 0.002812, p-value = 0.4078
    ## alternative hypothesis: two-sided

which confirms a good approximation.

## References

<div id="refs" class="references">

<div id="ref-Grzelak2014">

Grzelak, Lech A., Jeroen Witteveen, Maria Suarez-Taboada, and Cornelis
W. Oosterlee. 2014. “The Stochastic Collocation Monte Carlo Sampler:
Highly Efficient Sampling from ’Expensive’ Distributions.” *SSRN
Electronic Journal*. Elsevier BV.
<https://doi.org/10.2139/ssrn.2529691>.

</div>

</div>
