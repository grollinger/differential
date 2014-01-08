package rk

import (
	"errors"
	. "github.com/rollingthunder/differential/ode"
	"github.com/rollingthunder/differential/util"
	"math"
)

type RKMethod int

type rk struct {
	IntegratorInfo
	method           RKMethod
	firstStageAsLast bool
	b, c, e          []float64
	a                [][]float64
}

//-- performs Runge-Kutta integration
func (r *rk) Integrate(t, tEnd float64, yT []float64, c *Config) (stat Statistics, err error) {
	// set default parameters if necessary
	if c.MaxStepSize <= 0.0 {
		c.MaxStepSize = tEnd - t
	}
	if c.MinStepSize <= 0.0 {
		c.MinStepSize = 1e-10
	}
	if c.MaxStepCount == 0 {
		c.MaxStepCount = 1000000
	}
	if c.AbsoluteTolerance <= 0.0 {
		c.AbsoluteTolerance = 1e-4
	}
	if c.RelativeTolerance <= 0.0 {
		c.RelativeTolerance = c.AbsoluteTolerance
	}

	if r.a == nil || r.b == nil || r.c == nil {
		err = errors.New("RK Method coefficients not initialized")
		return
	}

	// local variables
	var n uint = uint(len(yT))

	// allocate temp matrices
	fcnValue := make([]float64, n)
	yCurrent := make([]float64, n)
	yError := make([]float64, n)
	ks := util.MakeRectangular(r.Stages, uint(n))

	c.Fcn(t, yT, fcnValue)
	stat.EvaluationCount = 1

	// compute initial step size if not set
	stepEstimate := c.InitialStepSize
	if stepEstimate <= 0.0 {
		stepEstimate = EstimateStepSize(t, yT, fcnValue, c, r.Order)
	}
	var stepNext float64
	// repeat until tend
	for t < tEnd && err == nil {
		// Set new step size
		stepNext = stepEstimate

		stat.StepCount++
		if t+stepNext > tEnd {
			stepNext = tEnd - t
		}

		// compute stages
		var stg, id, ic uint
		for stg = 1; stg < r.Stages; stg++ {
			tCurrent := t + stepNext*r.c[stg]

			for id = 0; id < n; id++ {
				yCurrent[id] = yT[id] + stepNext*r.a[stg][0]*fcnValue[id]
			}

			for ic = 1; ic < stg; ic++ {
				for id = 0; id < n; id++ {
					yCurrent[id] = yCurrent[id] + stepNext*r.a[stg][ic]*ks[ic][id]
				}
			}
			c.Fcn(tCurrent, yCurrent, ks[stg])
			stat.EvaluationCount++
		}

		// compute error estimate:
		for id = 0; id < n; id++ {
			yError[id] = stepNext * r.e[0] * fcnValue[id]
		}

		for stg = 1; stg < r.Stages; stg++ {
			for id = 0; id < n; id++ {
				yError[id] = yError[id] + stepNext*r.e[stg]*ks[stg][id]
			}
		}

		// compute error quotient
		relativeError := 0.0
		for id = 0; id < n; id++ {
			currentTolerance := c.AbsoluteTolerance + c.RelativeTolerance*math.Abs(yT[id])
			relativeError = relativeError + math.Pow(yError[id]/currentTolerance, 2.0)
		}
		relativeError = math.Sqrt(relativeError / float64(n))

		// new stepsize estimate
		stepEstimate = 0.9 * math.Exp(-math.Log(1.0e-8+relativeError)/float64(r.Order))
		stepEstimate = stepNext * math.Max(0.2, math.Min(stepEstimate, 2.0)) // safety interval

		// reject step
		if relativeError > 1.0 {
			stat.RejectedCount++

			// report failure, step size too small
			if stepEstimate < c.MinStepSize {
				err = errors.New("stepsize too small")
				break
			}
		} else {
			// accept step and compute new solution
			t += stepNext
			for id = 0; id < n; id++ {
				yT[id] = yT[id] + stepNext*r.b[0]*fcnValue[id]
			}
			for stg = 1; stg < r.Stages; stg++ {
				for id = 0; id < n; id++ {
					yT[id] = yT[id] + stepNext*r.b[stg]*ks[stg][id]
				}
			}

			// cancel after first step
			if c.OneStepOnly {
				break
			} else {
				if r.firstStageAsLast {
					copy(fcnValue, ks[r.Stages-1])
				} else {
					c.Fcn(t, yT, fcnValue)
					stat.EvaluationCount++
				}
			}
		}
		// failure, too many steps
		if stat.StepCount > c.MaxStepCount {
			err = errors.New("maximum step count exceeded")
			break
		}
	}

	stat.CurrentTime = t
	stat.LastStepSize = stepNext
	stat.NextStepSize = stepEstimate

	return

}
