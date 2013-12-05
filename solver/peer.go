package solver

import (
	"errors"
	"math"
)

type peer struct {
	IntegratorInfo
	method       PeerMethod
	icmin, icmax uint

	sigx, ema, emaoh float64
	c, e             []float64
	b, a0, cv, pv    [][]float64
}

func (p *peer) Integrate(t, tEnd float64, yT []float64, c Config) (stat Statistics, err error) {
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

	// local variables
	n := uint(len(yT))
	var old, newIdx int
	var tc, t0, rer, rc, sigmn, sig, sighs float64

	// allocate temp matrices
	yf := make([]float64, n)
	pa := makeSquare(p.stages)

	yy := [][][]float64{makeRectangular(p.stages, n), makeRectangular(p.stages, n)}
	ff := [][][]float64{makeRectangular(p.stages, n), makeRectangular(p.stages, n)}

	sigmn = 0.2

	// switch-indices old+newIdx=1
	old = 0
	newIdx = 1

	// startup with DOPRI
	dopri, err := NewRK(DoPri5)
	if err != nil {
		err = errors.New("error during startup: " + err.Error())
		return
	}

	copy(yy[old][p.icmin], yT)

	c.Fcn(t, yT, ff[old][p.icmin])
	stat.EvaluationCount = 1

	// guess initial step size if unspecified
	stepEstimate := c.InitialStepSize
	if stepEstimate <= 0.0 {
		stepEstimate = estimateStepSize(t, yT, ff[old][p.icmin], &c, p.order)
	}

	copy(yy[old][p.icmax], yT)
	tc = t

	//  higher accuracy for starting proc
	rkConfig := Config{
		InitialStepSize:   stepEstimate,
		RelativeTolerance: math.Max(1e-1*c.RelativeTolerance, 1e-14),
		AbsoluteTolerance: math.Max(1e-1*c.AbsoluteTolerance, 1e-14),
		OneStepOnly:       true,
		Fcn:               c.Fcn,
	}

	rkStat, err := dopri.Integrate(tc, t+stepEstimate, yy[old][p.icmax], rkConfig)
	if err != nil {
		err = errors.New("error during startup: " + err.Error())
		return
	}
	tc = rkStat.CurrentTime

	stat.EvaluationCount += rkStat.EvaluationCount
	c.Fcn(tc, yy[old][p.icmax], ff[old][p.icmax])
	stat.EvaluationCount++

	// adjusted step size, relative to interval [0,1]:
	stepPrevious := (tc - t) / (p.c[p.icmax] - p.c[p.icmin])
	t0 = t - stepPrevious*p.c[p.icmin] // corresponds to node pc=0

	// startup procedure
	var stg uint
	for stg = 0; stg < p.stages; stg++ {
		if stg != p.icmin && stg != p.icmax {
			copy(yy[old][stg], yT)
			rkConfig.InitialStepSize = stepPrevious * (p.c[stg] - p.c[p.icmin])
			tStage := t0 + stepPrevious*p.c[stg]
			rkStat, err = dopri.Integrate(t, tStage, yy[old][stg], rkConfig)
			if err != nil {
				err = errors.New("error during startup: " + err.Error())
				return
			}
			c.Fcn(tStage, yy[old][stg], ff[old][stg])
			stat.EvaluationCount += rkStat.EvaluationCount + 1
		}
	}

	// end of start interval
	t = t0 + stepPrevious
	stepEstimate = stepPrevious   // continue with stepsize stepPrevious
	copy(yT, yy[old][p.stages-1]) // solution at T

	var stepNext float64
	// repeat until tend
	for t < (tEnd - c.AbsoluteTolerance) {
		if t+stepEstimate > tEnd {
			stepEstimate = tEnd - t
		}
		stepNext = stepEstimate

		stat.StepCount++
		sig = stepNext / stepPrevious

		// COMPUTE COEFFS -> "Co" Prefix
		// stepPrevious*A row-wise
		// Loops: CoStages( CoA0, CoA1 )
		var stg, ic, id uint
		/*@; BEGIN(CoStages=Nest) @*/
		for stg = 0; stg < p.stages; stg++ {
			sighs = stepPrevious
			/*@; BEGIN(CoA0=Nest) @*/
			for ic = 0; ic < p.stages; ic++ {
				pa[stg][ic] = stepPrevious * p.a0[stg][ic]
			}
			/*@; BEGIN(CoA1=Nest) @*/
			for ic = 0; ic < p.stages; ic++ {
				sighs = sighs * sig
				for id = 0; id < p.stages; id++ {
					pa[stg][id] = pa[stg][id] + p.cv[stg][ic]*sighs*p.pv[ic][id]
				}
			}
		}

		// STAGE SOLUTIONS -> "St" Prefix
		// Loops: StB
		for id = 0; id < n; id++ {
			for stg = 0; stg < p.stages; stg++ {
				yy[newIdx][stg][id] = 0.0
				/*@; BEGIN(StB=Nest) @*/
				for ic = 0; ic < p.stages; ic++ {
					yy[newIdx][stg][id] = yy[newIdx][stg][id] + p.b[stg][ic]*yy[old][ic][id]
				}
			}
		}

		for id = 0; id < n; id++ {
			for stg = 0; stg < p.stages; stg++ {
				for ic = 0; ic < p.stages; ic++ {
					yy[newIdx][stg][id] = yy[newIdx][stg][id] + pa[stg][ic]*ff[old][ic][id]
				}
			}
		}

		// FUNCTION EVALUATIONS 
		// Fn=fcn(Yn)
		// Candidate for Parallelisation
		for stg = 0; stg < p.stages; stg++ {
			c.Fcn(t+stepNext*p.c[stg], yy[newIdx][stg], ff[newIdx][stg])
		}

		stat.EvaluationCount += p.stages

		// compute error estimate with F(newIdx):
		for id = 0; id < n; id++ {
			rc = 0.0
			for stg = 0; stg < p.stages; stg++ {
				rc = rc + p.e[stg]*ff[newIdx][stg][id]
			}
			yf[id] = math.Pow(rc/(c.AbsoluteTolerance+c.RelativeTolerance*math.Abs(yT[id])), 2.0)
		}

		// compute error quotient/20070803
		// step ratio from error model ((1+a)^p-a^p)/est+a^p)^(1/p)-a, p=order/2:
		rer = 0.0
		for id = 0; id < n; id++ {
			rer = rer + yf[id]
		}

		rer = stepNext*math.Sqrt(rer/float64(n)) + 1e-8
		sighs = math.Pow(math.Pow(sig, 2.0)+p.ema, float64(p.order)/2.0) - p.emaoh
		stepEstimate = math.Pow(sighs/rer+p.emaoh, 2.0/float64(p.order)) - p.ema
		stepEstimate = stepPrevious * math.Max(sigmn, math.Min(0.95*math.Sqrt(stepEstimate), p.sigx)) // safety interval

		// reject step
		if rer > 1.0 {
			// decrease minimal stepsize ratio
			sigmn = sigmn * 0.2
			stat.RejectedCount++

			// report failure
			if stepEstimate < c.MinStepSize {
				err = errors.New("step size too small")
				break
			}
		} else {
			// accept step, swap YY & FF
			sigmn = 0.2
			old = newIdx
			newIdx = 1 - old
			copy(yT, yy[old][p.stages-1])
			t = t + stepNext
			stepPrevious = stepNext
		}

		// failure, too many steps
		if stat.StepCount > c.MaxStepCount {
			err = errors.New("maximum step count exceeded")
			break
		}
	}

	// set return code and time
	stat.CurrentTime = t
	stat.LastStepSize = stepNext
	stat.NextStepSize = stepEstimate
	return
}
