// Copyright ©2017 The Gonum Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mat

import (
	"testing"

	"golang.org/x/exp/rand"

	"gonum.org/v1/gonum/floats"
)

func TestGSVD(t *testing.T) {
	const tol = 1e-10
	rnd := rand.New(rand.NewSource(1))
	for _, test := range []struct {
		m, p, n int
	}{
		{5, 3, 5},
		{5, 3, 3},
		{3, 3, 5},
		{5, 5, 5},
		{5, 5, 3},
		{3, 5, 5},
		{150, 150, 150},
		{200, 150, 150},
		{150, 150, 200},
		{150, 200, 150},
		{200, 200, 150},
		{150, 200, 200},
	} {
		m := test.m
		p := test.p
		n := test.n
		for trial := 0; trial < 10; trial++ {
			a := NewDense(m, n, nil)
			for i := range a.mat.Data {
				a.mat.Data[i] = rnd.NormFloat64()
			}
			aCopy := DenseCopyOf(a)

			b := NewDense(p, n, nil)
			for i := range b.mat.Data {
				b.mat.Data[i] = rnd.NormFloat64()
			}
			bCopy := DenseCopyOf(b)

			// Test Full decomposition.
			var gsvd GSVD
			ok := gsvd.Factorize(a, b)
			if !ok {
				t.Errorf("GSVD factorization failed")
			}
			if !Equal(a, aCopy) {
				t.Errorf("A changed during call to GSVD.Factorize with full decomposition")
			}
			if !Equal(b, bCopy) {
				t.Errorf("B changed during call to GSVD.Factorize with full decomposition")
			}
			c, s, sigma1, sigma2, zeroR, u, v, q := extractGSVD(&gsvd)
			var ansU, ansV, d1R, d2R Dense
			ansU.Product(u.T(), a, q)
			ansV.Product(v.T(), b, q)
			d1R.Mul(sigma1, zeroR)
			d2R.Mul(sigma2, zeroR)
			if !EqualApprox(&ansU, &d1R, tol) {
				t.Errorf("Answer mismatch with full decomposition\nU^T * A * Q:\n% 0.2f\nΣ₁ * [ 0 R ]:\n% 0.2f",
					Formatted(&ansU), Formatted(&d1R))
			}
			if !EqualApprox(&ansV, &d2R, tol) {
				t.Errorf("Answer mismatch with full decomposition\nV^T * B  *Q:\n% 0.2f\nΣ₂ * [ 0 R ]:\n% 0.2f",
					Formatted(&d2R), Formatted(&ansV))
			}

			// Check C^2 + S^2 = I.
			for i := range c {
				d := c[i]*c[i] + s[i]*s[i]
				if !floats.EqualWithinAbsOrRel(d, 1, 1e-14, 1e-14) {
					t.Errorf("c_%d^2 + s_%d^2 != 1: got: %v", i, i, d)
				}
			}

			// Test None decomposition.
			gsvd.NoU = true
			gsvd.NoV = true
			gsvd.NoQ = true
			ok = gsvd.Factorize(a, b)
			if !ok {
				t.Errorf("GSVD factorization failed")
			}
			if !Equal(a, aCopy) {
				t.Errorf("A changed during call to GSVD with no vectors")
			}
			if !Equal(b, bCopy) {
				t.Errorf("B changed during call to GSVD with no vectors")
			}
			cNone := gsvd.ValuesA(nil)
			if !floats.EqualApprox(c, cNone, tol) {
				t.Errorf("Singular value mismatch between full and no decomposition")
			}
			sNone := gsvd.ValuesB(nil)
			if !floats.EqualApprox(s, sNone, tol) {
				t.Errorf("Singular value mismatch between full and no decomposition")
			}
		}
	}
}

func extractGSVD(gsvd *GSVD) (c, s []float64, s1, s2, zR, u, v, q *Dense) {
	s1 = gsvd.SigmaATo(nil)
	s2 = gsvd.SigmaBTo(nil)
	zR = gsvd.ZeroRTo(nil)
	u = gsvd.UTo(nil)
	v = gsvd.VTo(nil)
	q = gsvd.QTo(nil)
	c = gsvd.ValuesA(nil)
	s = gsvd.ValuesB(nil)
	return c, s, s1, s2, zR, u, v, q
}
