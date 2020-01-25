// Copyright ©2014 The Gonum Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package distuv

import (
	"math"

	"golang.org/x/exp/rand"
)

// #cgo LDFLAGS: -lgsl -lgslcblas
// #include "./libs/libstable/stable/src/mcculloch.c"
// #include "./libs/libstable/stable/src/mcculloch.h"
// #include "./libs/libstable/stable/src/methods.c"
// #include "./libs/libstable/stable/src/methods.h"
// #include "./libs/libstable/stable/src/stable_cdf.c"
// #include "./libs/libstable/stable/src/stable_common.c"
// #include "./libs/libstable/stable/src/stable_dist.c"
// #include "./libs/libstable/stable/src/stable_fit.c"
// #include "./libs/libstable/stable/src/stable_fit.h"
// #include "./libs/libstable/stable/src/stable_integration.c"
// #include "./libs/libstable/stable/src/stable_integration.h"
// #include "./libs/libstable/stable/src/stable_koutrouvelis.c"
// #include "./libs/libstable/stable/src/stable_pdf.c"
// #include "./libs/libstable/stable/src/stable_q.c"
// #include "./libs/libstable/stable/src/stable_rnd.c"
// #include "./libs/libstable/stable/src/stable.h"
import "C"

// StandardStable returns an instantiation of the stable distribution with Mu = 0 and Sigma = 1.
func StandardStable(alpha float64, beta float64) Stable {
	return Stable{Alpha: alpha, Beta: beta, Mu: 0, Sigma: 1}
}

// Stable respresents a Levy alpha-stable distribution (https://en.wikipedia.org/wiki/Stable_distribution).
type Stable struct {
	Alpha float64
	Beta  float64
	Mu    float64 // Location param of the stable distribution (Mean for alpha > 1)
	Sigma float64 // Scale param of the stable distribution (StdDev/sqrt2 if alpha == 2; Normal distribution)
	Src   rand.Source
}

// CDF computes the value of the cumulative density function at x.
//
// TODO:
func (s Stable) CDF(x float64) float64 {
	return 0.0
}

// ConjugateUpdate updates the parameters of the distribution from the sufficient
// statistics of a set of samples. The sufficient statistics, suffStat, have been
// observed with nSamples observations. The prior values of the distribution are those
// currently in the distribution, and have been observed with priorStrength samples.
//
// TODO:
func (s *Stable) ConjugateUpdate(suffStat []float64, nSamples float64, priorStrength []float64) {
}

// Entropy returns the differential entropy of the distribution.
//
// TODO:
func (s Stable) Entropy() float64 {
	return 0.0
}

// ExKurtosis returns the excess kurtosis of the distribution.
func (s Stable) ExKurtosis() float64 {
	if s.Alpha == 2 {
		return 0
	}
	return math.NaN()
}

// Fit sets the parameters of the probability distribution from the
// data samples x with relative weights w. If weights is nil, then all the weights
// are 1. If weights is not nil, then the len(weights) must equal len(samples).
//
// Fit to data using: http://auapps.american.edu/jpnolan/www/stable/mle.pdf
// Adapted from libstable/matlab/stable_fit_mle2dC.m
// NOTE: Be sure to install gsl (e.g. brew install gsl)
// For references, see http://www.lpi.tel.uva.es/stable and https://golang.org/cmd/cgo/
func (s *Stable) Fit(samples, weights []float64) {
	// NOTE: Ignoring weights for now ...
	// Safely build the C.double version of samples slice
	n := len(samples)
	cs := make([]C.double, n)
	for i, sample := range samples {
		cs[i] = C.double(sample)
	}

	// Create initial stable dist from existing s params
	dist := C.stable_create(C.double(s.Alpha), C.double(s.Beta), C.double(s.Sigma), C.double(s.Mu), 1) // Parametrization is 1 in Chambers, J.M., Mallows, C.L. and Stuck, B.W.

	// Set config vals
	C.stable_set_THREADS(0)
	C.stable_set_relTOL(1e-8)
	C.stable_set_absTOL(1e-8)

	// Fit iterative to samples
	if status := C.stable_fit(dist, &(cs[0]), C.uint(n)); int(status) > 0 {
		panic("stable fit: an error occured on iterative procedure")
	}

	// Set the fitted values for params
	s.Alpha = float64(dist.alpha)
	s.Beta = float64(dist.beta)
	s.Mu = float64(dist.mu_1)
	s.Sigma = float64(dist.sigma)
}

// LogProb computes the natural logarithm of the value of the probability density function at x.
//
// TODO:
func (s Stable) LogProb(x float64) float64 {
	return 0
}

// Mean returns the mean of the probability distribution.
func (s Stable) Mean() float64 {
	if s.Alpha > 1 {
		return s.Mu
	}
	return math.NaN()
}

// Median returns the median of the stable distribution.
func (s Stable) Median() float64 {
	if s.Beta == 0 {
		return s.Mu
	}
	return math.NaN()
}

// Mode returns the mode of the stable distribution.
func (s Stable) Mode() float64 {
	if s.Beta == 0 {
		return s.Mu
	}
	return math.NaN()
}

// NumParameters returns the number of parameters in the distribution.
func (Stable) NumParameters() int {
	return 4
}

// NumSuffStat returns the number of sufficient statistics for the distribution.
//
// TODO:
func (Stable) NumSuffStat() int {
	return 4
}

// Prob computes the value of the probability density function at x.
//
// TODO:
func (s Stable) Prob(x float64) float64 {
	return 0
}

// Quantile returns the inverse of the cumulative probability distribution.
//
// Adapted from libstable/matlab/stable_invC.m
// NOTE: Be sure to install gsl (e.g. brew install gsl)
// For references, see http://www.lpi.tel.uva.es/stable and https://golang.org/cmd/cgo/
func (s Stable) Quantile(p float64) float64 {
	if p < 0 || p > 1 {
		panic(badPercentile)
	}

	// Create initial stable dist from existing s params
	dist := C.stable_create(C.double(s.Alpha), C.double(s.Beta), C.double(s.Sigma), C.double(s.Mu), 1) // Parametrization is 1 in Chambers, J.M., Mallows, C.L. and Stuck, B.W.

	// Set config vals
	C.stable_set_THREADS(0)
	C.stable_set_relTOL(1e-12)
	C.stable_set_absTOL(1e-16)
	C.stable_set_INV_MAXITER(50)

	// Determine and return CDF^-1 value
	return float64(C.stable_q_point(dist, C.double(p), nil))
}

// Rand returns a random sample drawn from the distribution.
func (s Stable) Rand() float64 {
	// TODO: Panic for invalid values of alpha, beta, sigma

	// Generate using
	//  Chambers, J.M., Mallows, C.L. and Stuck, B.W.. "A Method for Simulating
	//  Stable Random Variables"
	//  https://doi.org/10.1080%2F01621459.1976.10480344
	//  use this reference: https://www.tandfonline.com/doi/abs/10.1080/01621459.1976.10480344

	u := Uniform{Min: -math.Pi / 2, Max: math.Pi / 2, Src: s.Src}.Rand()
	w := Exponential{Rate: 1, Src: s.Src}.Rand()
	zeta := -s.Beta * math.Tan(math.Pi*s.Alpha/2)
	if s.Alpha != 1 {
		xi := (1 / s.Alpha) * math.Atan(-zeta)
		stndrnd := math.Pow(1+math.Pow(zeta, 2), 1/(2*s.Alpha)) * (math.Sin(s.Alpha*(u+xi)) / math.Pow(math.Cos(u), 1/s.Alpha)) * math.Pow(math.Cos(u-s.Alpha*(u+xi))/w, (1-s.Alpha)/s.Alpha)
		return s.Sigma*stndrnd + s.Mu
	} else {
		xi := math.Pi / 2
		stndrnd := (1 / xi) * ((xi+s.Beta*u)*math.Tan(u) - s.Beta*math.Log(w*math.Cos(u)/(xi+s.Beta*u)))
		return s.Sigma*stndrnd + s.Mu + (2/math.Pi)*s.Beta*s.Sigma*math.Log(s.Sigma)
	}
}

// Score returns the score function with respect to the parameters of the
// distribution at the input location x. The score function is the derivative
// of the log-likelihood at x with respect to the parameters
//  (∂/∂θ) log(p(x;θ))
// If deriv is non-nil, len(deriv) must equal the number of parameters otherwise
// Score will panic, and the derivative is stored in-place into deriv. If deriv
// is nil a new slice will be allocated and returned.
//
// The order is [∂LogProb / ∂Mu, ∂LogProb / ∂Sigma].
//
// For more information, see https://en.wikipedia.org/wiki/Score_%28statistics%29.
//
// TODO:
func (s Stable) Score(deriv []float64, x float64) []float64 {
	return deriv
}

// ScoreInput returns the score function with respect to the input of the
// distribution at the input location specified by x. The score function is the
// derivative of the log-likelihood
//  (d/dx) log(p(x)) .
func (s Stable) ScoreInput(x float64) float64 {
	return 0
}

// Skewness returns the skewness of the distribution.
func (s Stable) Skewness() float64 {
	if s.Alpha == 2 {
		return 0
	}
	return math.NaN()
}

// StdDev returns the standard deviation of the probability distribution.
func (s Stable) StdDev() float64 {
	if s.Alpha == 2 {
		return math.Sqrt2 * s.Sigma
	}
	return math.NaN()
}

// SuffStat computes the sufficient statistics of a set of samples to update
// the distribution. The sufficient statistics are stored in place, and the
// effective number of samples are returned.
//
// TODO:
//
// If weights is nil, the weights are assumed to be 1, otherwise panics if
// len(samples) != len(weights). Panics if len(suffStat) != NumSuffStat().
func (Stable) SuffStat(suffStat, samples, weights []float64) (nSamples float64) {
	nSamples = 0
	return nSamples
}

// Survival returns the survival function (complementary CDF) at x.
//
// TODO:
func (s Stable) Survival(x float64) float64 {
	return 0
}

// SetParameters modifies the parameters of the distribution.
func (s *Stable) SetParameters(p []Parameter) {
	if len(p) != s.NumParameters() {
		panic("stable: incorrect number of parameters to set")
	}
	if p[0].Name != "Alpha" {
		panic("stable: " + panicNameMismatch)
	}
	if p[1].Name != "Beta" {
		panic("stable: " + panicNameMismatch)
	}
	if p[2].Name != "Mu" {
		panic("stable: " + panicNameMismatch)
	}
	if p[3].Name != "Sigma" {
		panic("stable: " + panicNameMismatch)
	}
	s.Alpha = p[0].Value
	s.Beta = p[1].Value
	s.Mu = p[2].Value
	s.Sigma = p[3].Value
}

// Variance returns the variance of the probability distribution.
func (s Stable) Variance() float64 {
	if s.Alpha == 2 {
		return 2 * math.Pow(s.Sigma, 2)
	}
	return math.NaN()
}

// parameters returns the parameters of the distribution.
func (s Stable) parameters(p []Parameter) []Parameter {
	nParam := s.NumParameters()
	if p == nil {
		p = make([]Parameter, nParam)
	} else if len(p) != nParam {
		panic("stable: improper parameter length")
	}
	p[0].Name = "Alpha"
	p[0].Value = s.Alpha
	p[1].Name = "Beta"
	p[1].Value = s.Beta
	p[2].Name = "Mu"
	p[2].Value = s.Mu
	p[3].Name = "Sigma"
	p[3].Value = s.Sigma
	return p
}
