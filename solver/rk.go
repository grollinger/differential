package solver

import (
	"errors"
	"math"
)

type RKMethod int

type rk struct {
	parameters
	method           RKMethod
	firstStageAsLast bool
	b, c, e          []float64
	a                [][]float64
}

type Statistics struct {
	StepCount, RejectedCount, EvaluationCount int
	LastStepSize, CurrentTime                 float64
}

func (r *rk) SetParameters(atol, rtol float64) {
	r.atol, r.rtol = atol, rtol
}

func (r *rk) IntegrateStep(stepSize, t float64, y []float64, fcn Function) {

}

//-- performs Runge-Kutta integration
func (r *rk) Integrate(initialStepSize, t, tEnd float64, y []float64, fcn Function) (stat Statistics, err error) {
	// local variables
	var n uint = uint(len(y))

	var tc, rer, rc, hn float64
	hs := initialStepSize

	var o1step bool

	// allocate temp matrices
	fv1 := make([]float64, n)
	yc := make([]float64, n)
	ye := make([]float64, n)
	ks := makeRectangular(r.stages, uint(n))

	r.hmax = tEnd - t

	// rksetcoeff not called
	if r.a == nil {
		err = errors.New("rk: method coefficients not initialized")
		return
	}

	fcn(t, y, fv1)
	stat.EvaluationCount = 1

	// compute initial stepsize
	if hs <= 0.0 {
		hs = estimateStepsize(t, y, fv1, fcn, r.parameters)
	}

	// repeat until tend
	for t < tEnd && err == nil {
		//printf("Step %d : y = %1.20e \t %1.20e\n", rksteps, y[0], hs);

		stat.StepCount++
		if t+hs > tEnd {
			hs = tEnd - t
		}

		// compute stages
		var stg, id, ic uint
		for stg = 1; stg < r.stages; stg++ {
			tc = t + hs*r.c[stg]

			for id = 0; id < n; id++ {
				yc[id] = y[id] + hs*r.a[stg][0]*fv1[id]
			}

			for ic = 1; ic < stg; ic++ {
				for id = 0; id < n; id++ {
					yc[id] = yc[id] + hs*r.a[stg][ic]*ks[ic][id]
				}
			}
			fcn(tc, yc, ks[stg])
			stat.EvaluationCount++
		}

		// compute error estimate:
		for id = 0; id < n; id++ {
			ye[id] = hs * r.e[0] * fv1[id]
		}

		for stg = 1; stg < r.stages; stg++ {
			for id = 0; id < n; id++ {
				ye[id] = ye[id] + hs*r.e[stg]*ks[stg][id]
			}
		}

		// compute error quotient
		rer = 0.0
		for id = 0; id < n; id++ {
			rc = r.atol + r.rtol*math.Abs(y[id])
			rer = rer + math.Pow(ye[id]/rc, 2.0)
		}
		rer = math.Sqrt(rer / float64(n))

		// new stepsize estimate
		hn = 0.9 * math.Exp(-math.Log(1.0e-8+rer)/float64(r.order))
		hn = hs * math.Max(0.2, math.Min(hn, 2.0)) // safety interval

		// reject step
		if rer > 1.0 {
			stat.RejectedCount++

			// report failure, stepsize too small
			if hn < r.hmin {
				err = errors.New("rk: stepsize too small")
				return
			}
		} else {
			// accept step and compute new solution
			t += hs
			for id = 0; id < n; id++ {
				y[id] = y[id] + hs*r.b[0]*fv1[id]
			}
			for stg = 1; stg < r.stages; stg++ {
				for id = 0; id < n; id++ {
					y[id] = y[id] + hs*r.b[stg]*ks[stg][id]
				}
			}

			// cancel after first step
			if o1step {
				//TODO remove
				return
			} else {
				if r.firstStageAsLast {
					copy(fv1, ks[r.stages-1])
				} else {
					fcn(t, y, fv1)
					stat.EvaluationCount++
				}
			}
		}
		// failure, too many steps
		if stat.StepCount > r.maxsteps {
			err = errors.New("rk: too many steps")
			return
		}

		// Set new stepsize
		hs = hn

	}

	// set return code and time

	stat.CurrentTime = t
	stat.LastStepSize = hs

	return

}
