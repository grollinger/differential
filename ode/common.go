package ode

import "math"

func EstimateStepSize(t float64, yT, fcnValue []float64, c *Config, order uint) float64 {
	n := len(yT)
	var h, h1, der2, der12 float64

	// allocate temp arrays
	y2, f2 := make([]float64, n), make([]float64, n)

	// calculate temp step size
	dnf, dny := 0.0, 0.0
	for id := 0; id < n; id++ {
		rc := c.AbsoluteTolerance + c.RelativeTolerance*math.Abs(yT[id])
		dnf = dnf + math.Pow(fcnValue[id]/rc, 2)
		dny = dny + math.Pow(yT[id]/rc, 2)
	}

	if math.Min(dnf, dny) < 1e-10 {
		h = 1.e-6
	} else {
		h = 1.e-2 * math.Sqrt(dny/dnf)
	}
	h = math.Min(h, c.MaxStepSize)

	// explicit Euler step
	for id := 0; id < n; id++ {
		y2[id] = yT[id] + h*fcnValue[id]
	}
	c.Fcn(t+h, y2, f2)

	der2 = 0.0
	for id := 0; id < n; id++ {
		rc := c.AbsoluteTolerance + c.RelativeTolerance*math.Abs(yT[id])
		der2 = der2 + math.Pow((f2[id]-fcnValue[id])/rc, 2)
	}

	//estimate for second derivative
	der2 = math.Sqrt(der2) / h
	der12 = math.Max(der2, math.Sqrt(dnf))

	// calculate initial stepsize
	if der12 <= 1.e-15 {
		h1 = math.Max(1.e-6, h*1.e-3)
	} else {
		h1 = math.Pow(1.e-2/der12, 1.0/float64(order))
	}
	return math.Min(1e2*h, math.Min(h1, c.MaxStepSize))
}

//-- solves Vandermonde systems PM_new*V=PM
func VanderMonde(pc []float64, pm [][]float64) {
	var i, j, k int
	n := len(pc)

	for k = 0; k < n-1; k++ {
		for j = n - 1; j >= k+1; j-- {
			for i = 0; i < n; i++ {
				pm[i][j] = pm[i][j] - pm[i][j-1]*pc[k]
			}
		}
	}

	for k = n - 2; k >= 0; k-- {
		for j = k + 1; j < n; j++ {
			for i = 0; i < n; i++ {
				pm[i][j] = pm[i][j] / (pc[j] - pc[j-k-1])
			}
		}

		for j = k; j < n-1; j++ {
			for i = 0; i < n; i++ {
				pm[i][j] = pm[i][j] - pm[i][j+1]
			}
		}
	}
}
