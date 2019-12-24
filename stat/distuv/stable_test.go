// Copyright ©2017 The Gonum Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package distuv

import (
	"log"
	"math"
	"testing"

	"golang.org/x/exp/rand"
)

func TestStableProb(t *testing.T) {
	const tol = 1e-10
}

func TestStableCDF(t *testing.T) {
	const tol = 1e-10
}

func TestStable(t *testing.T) {
}

func TestStableFit(t *testing.T) {
	//   1. Init a stable dist with certain param vals (alpha, beta, mu, sigma)
	//   2. Generate n random values from dist for samples
	//   3. Call Fit function wrt prev generated sample values
	//   4. Check fitted params are close in value to orig dist params
	// Pick init dist
	alpha := 1.4
	beta := 0.1
	mu := 0.02
	sigma := 0.05
	src := rand.New(rand.NewSource(1000))
	smplDist := Stable{
		Alpha: alpha,
		Beta:  beta,
		Mu:    mu,
		Sigma: sigma,
		Src:   src,
	}

	n := 5000
	samples := make([]float64, n)
	for i := 0; i < n; i++ {
		samples[i] = smplDist.Rand()
	}

	// Init fitting dist
	dist := Stable{
		Alpha: 1.0,
		Beta:  0.0,
		Mu:    0.0,
		Sigma: 1.0,
		Src:   src,
	}
	log.Println(dist)
	dist.Fit(samples, []float64{})
	log.Println(dist)

	// TODO: Check param values are within tolerance range
}

func TestStableQuantile(t *testing.T) {
	// Pick init dist
	alpha := 1.4
	beta := 0.1
	mu := 0.02
	sigma := 0.05
	src := rand.New(rand.NewSource(1000))
	dist := Stable{
		Alpha: alpha,
		Beta:  beta,
		Mu:    mu,
		Sigma: sigma,
		Src:   src,
	}
	z := dist.Quantile(0.99)
	want := 0.53237
	tol := 1e-6

	if math.Abs(want-z) > tol {
		t.Errorf("test-quantile: got=%e. want=%e", z, want)
	}
}
