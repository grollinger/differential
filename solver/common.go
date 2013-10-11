package solver

import "math"

type Function func(float64, []float64, []float64)

type parameters struct {
	stages, order uint
	//misc params
	maxsteps               int
	hmin, hmax, atol, rtol float64
}

func makeSquare(n uint) [][]float64 {
	return makeRectangular(n, n)
}

func makeRectangular(rows, cols uint) (rect [][]float64) {
	arr := make([]float64, rows*cols)
	rect = make([][]float64, rows)
	for i := range rect {
		rect[i] = arr[:cols]
		arr = arr[cols:]
	}
	return
}

func estimateStepsize(t float64, y, fv1 []float64, fcn Function, p parameters) float64 {
	n := len(y)
	var h, h1, rc, der2, der12 float64

	// allocate temp arrays
	y2, f2 := make([]float64, n), make([]float64, n)

	// calculate temp stepsize
	dnf, dny := 0.0, 0.0
	for id := 0; id < n; id++ {
		rc = p.atol + p.rtol*math.Abs(y[id])
		dnf = dnf + math.Pow(fv1[id]/rc, 2)
		dny = dny + math.Pow(y[id]/rc, 2)
	}

	if math.Min(dnf, dny) < 1e-10 {
		h = 1.e-6
	} else {
		h = 1.e-2 * math.Sqrt(dny/dnf)
	}
	h = math.Min(h, p.hmax)

	// explicit Euler step
	for id := 0; id < n; id++ {
		y2[id] = y[id] + h*fv1[id]
	}
	fcn(t+h, y2, f2)

	der2 = 0.0
	for id := 0; id < n; id++ {
		rc = p.atol + p.rtol*math.Abs(y[id])
		der2 = der2 + math.Pow((f2[id]-fv1[id])/rc, 2)
	}

	//estimate for second derivative
	der2 = math.Sqrt(der2) / h
	der12 = math.Max(der2, math.Sqrt(dnf))

	// calculate initial stepsize
	if der12 <= 1.e-15 {
		h1 = math.Max(1.e-6, h*1.e-3)
	} else {
		h1 = math.Pow(1.e-2/der12, 1.0/float64(p.order))
	}
	return math.Min(1e2*h, math.Min(h1, p.hmax))
}
