package epp

import (
	"errors"
	. "github.com/rollingthunder/differential/ode"
	"github.com/rollingthunder/differential/ode/rk"
	"github.com/rollingthunder/differential/util"
	"math"
)

type peer struct {
	IntegratorInfo
	method PeerMethod
	peerCoefficients
}

type peerCoefficients struct {
	indexMinNode, indexMaxNode uint

	stepRatioMax float64

	// Error Model
	// p=order/2
	// ((1+a)^p-a^p)/est+a^p)^(1/p)-a
	errorModelA float64
	// a0 = a^p
	errorModelA0      float64
	errorModelWeights []float64

	c             []float64
	b, a0, cv, pv [][]float64
}

type integration struct {
	Config
	Statistics
	fOld, fNew, yOld, yNew, pa                              [][]float64
	yf                                                      []float64
	t0, stepRatioMin, stepRatio, stepEstimate, stepPrevious float64
	n                                                       uint
}

func (p *peer) Integrate(t, tEnd float64, yT []float64, cfg Config) (stat Statistics, err error) {
	in := p.setupIntegration(t, tEnd, yT, cfg)

	p.startupIntegration(&in)

	// end of start interval
	t = in.t0 + in.stepPrevious
	in.stepEstimate = in.stepPrevious // continue with stepsize stepPrevious

	var stepNext float64
	// repeat until tend
	for t < (tEnd - in.AbsoluteTolerance) {
		if t+in.stepEstimate > tEnd {
			in.stepEstimate = tEnd - t
		}
		stepNext = in.stepEstimate

		stat.StepCount++
		in.stepRatio = stepNext / in.stepPrevious

		p.computeCoefficients(&in)

		// STAGE SOLUTIONS -> "St" Prefix
		var stg, ic, id uint
		// Loops: StB
		for id = 0; id < in.n; id++ {
			for stg = 0; stg < p.Stages; stg++ {
				in.yNew[stg][id] = 0.0
				/*@; BEGIN(StB=Nest) @*/
				for ic = 0; ic < p.Stages; ic++ {
					in.yNew[stg][id] += p.b[stg][ic] * in.yOld[ic][id]
				}
			}
		}

		for id = 0; id < in.n; id++ {
			for stg = 0; stg < p.Stages; stg++ {
				for ic = 0; ic < p.Stages; ic++ {
					in.yNew[stg][id] += in.pa[stg][ic] * in.fOld[ic][id]
				}
			}
		}

		// FUNCTION EVALUATIONS
		// Fn=fcn(Yn)
		// Candidate for Parallelisation
		for stg = 0; stg < p.Stages; stg++ {
			in.Fcn(t+stepNext*p.c[stg], in.yNew[stg], in.fNew[stg])
		}

		stat.EvaluationCount += p.Stages

		// compute error estimate with fNew:
		for id = 0; id < in.n; id++ {
			rc := 0.0
			for stg = 0; stg < p.Stages; stg++ {
				rc += p.errorModelWeights[stg] * in.fNew[stg][id]
			}
			in.yf[id] = math.Pow(rc/(in.AbsoluteTolerance+in.RelativeTolerance*math.Abs(in.yOld[p.Stages-1][id])), 2.0)
		}

		// compute error quotient/20070803
		// step ratio from error model ((1+a)^p-a^p)/est+a^p)^(1/p)-a, p=order/2:
		errorRelative := 0.0
		for id = 0; id < in.n; id++ {
			errorRelative += in.yf[id]
		}

		errorEstimate := stepNext*math.Sqrt(errorRelative/float64(in.n)) + 1e-8
		errorModelDenom := math.Pow(math.Pow(in.stepRatio, 2.0)+p.errorModelA, float64(p.Order)/2.0) - p.errorModelA0
		in.stepEstimate = math.Pow(errorModelDenom/errorEstimate+p.errorModelA0, 2.0/float64(p.Order)) - p.errorModelA
		in.stepEstimate = in.stepPrevious * math.Max(in.stepRatioMin, math.Min(0.95*math.Sqrt(in.stepEstimate), p.stepRatioMax)) // safety interval

		if errorEstimate > 1.0 {
			// reject step
			// decrease minimal stepsize ratio
			in.stepRatioMin = in.stepRatioMin * 0.2
			stat.RejectedCount++

			// report failure
			if in.stepEstimate < in.MinStepSize {
				err = errors.New("step size too small")
				break
			}
		} else {
			// accept step
			in.stepRatioMin = 0.2

			// swap Y & F
			swap := in.yOld
			in.yOld = in.yNew
			in.yNew = swap

			swap = in.fOld
			in.fOld = in.fNew
			in.fNew = swap

			t = t + stepNext
			in.stepPrevious = stepNext
		}

		// failure, too many steps
		if stat.StepCount > in.MaxStepCount {
			err = errors.New("maximum step count exceeded")
			break
		}
	}
	// output of last stage is the final output
	copy(yT, in.yOld[p.Stages-1])

	stat.CurrentTime = t
	stat.LastStepSize = in.stepPrevious
	stat.NextStepSize = in.stepEstimate
	return
}

func (p *peer) setupIntegration(t, tEnd float64, yT []float64, c Config) (i integration) {
	i.Config = c

	// set default parameters if necessary
	if i.MaxStepSize <= 0.0 {
		i.MaxStepSize = tEnd - t
	}
	if i.MinStepSize <= 0.0 {
		i.MinStepSize = 1e-10
	}
	if i.MaxStepCount == 0 {
		i.MaxStepCount = 1000000
	}
	if i.AbsoluteTolerance <= 0.0 {
		i.AbsoluteTolerance = 1e-4
	}
	if i.RelativeTolerance <= 0.0 {
		i.RelativeTolerance = i.AbsoluteTolerance
	}

	// local variables
	i.n = uint(len(yT))

	// allocate temp matrices
	i.yf = make([]float64, i.n)
	i.pa = util.MakeSquare(p.Stages)

	i.yNew = util.MakeRectangular(p.Stages, i.n)
	i.yOld = util.MakeRectangular(p.Stages, i.n)
	i.fNew = util.MakeRectangular(p.Stages, i.n)
	i.fOld = util.MakeRectangular(p.Stages, i.n)

	copy(i.yOld[p.indexMinNode], yT)

	i.stepRatioMin = 0.2
	i.t0 = t

	return
}

func (p *peer) startupIntegration(in *integration) {
	// startup with DOPRI
	dopri, err := rk.NewRK(rk.DoPri5)
	if err != nil {
		err = errors.New("error during startup: " + err.Error())
		return
	}

	in.Fcn(in.t0, in.yOld[p.indexMinNode], in.fOld[p.indexMinNode])
	in.EvaluationCount = 1

	// guess initial step size if unspecified
	in.stepEstimate = in.InitialStepSize
	if in.stepEstimate <= 0.0 {
		in.stepEstimate = EstimateStepSize(in.t0, in.yOld[p.indexMinNode], in.fOld[p.indexMinNode], &in.Config, p.Order)
	}

	copy(in.yOld[p.indexMaxNode], in.yOld[p.indexMinNode])
	tCurrent := in.t0

	//  higher accuracy for starting proc
	rkConfig := Config{
		InitialStepSize:   in.stepEstimate,
		RelativeTolerance: math.Max(1e-1*in.RelativeTolerance, 1e-14),
		AbsoluteTolerance: math.Max(1e-1*in.AbsoluteTolerance, 1e-14),
		OneStepOnly:       true,
		Fcn:               in.Fcn,
	}

	rkStat, err := dopri.Integrate(tCurrent, in.t0+in.stepEstimate, in.yOld[p.indexMaxNode], rkConfig)
	if err != nil {
		err = errors.New("error during startup: " + err.Error())
		return
	}
	tCurrent = rkStat.CurrentTime

	in.EvaluationCount += rkStat.EvaluationCount
	in.Fcn(tCurrent, in.yOld[p.indexMaxNode], in.fOld[p.indexMaxNode])
	in.EvaluationCount++

	// adjusted step size, relative to interval [0,1]:
	in.stepPrevious = (tCurrent - in.t0) / (p.c[p.indexMaxNode] - p.c[p.indexMinNode])
	origT0 := in.t0
	in.t0 -= in.stepPrevious * p.c[p.indexMinNode] // corresponds to node pc=0

	// startup procedure
	var stg uint
	for stg = 0; stg < p.Stages; stg++ {
		if stg != p.indexMinNode && stg != p.indexMaxNode {
			copy(in.yOld[stg], in.yOld[p.indexMinNode])
			rkConfig.InitialStepSize = in.stepPrevious * (p.c[stg] - p.c[p.indexMinNode])
			tStage := in.t0 + in.stepPrevious*p.c[stg]
			rkStat, err = dopri.Integrate(origT0, tStage, in.yOld[stg], rkConfig)
			if err != nil {
				err = errors.New("error during startup: " + err.Error())
				return
			}
			in.Fcn(tStage, in.yOld[stg], in.fOld[stg])
			in.EvaluationCount += rkStat.EvaluationCount + 1
		}
	}
}

func (p *peer) computeCoefficients(in *integration) {
	// COMPUTE COEFFS -> "Co" Prefix
	// stepPrevious*A row-wise
	// Loops: CoStages( CoA0, CoA1 )
	var stg, ic, id uint
	/*@; BEGIN(CoStages=Nest) @*/
	for stg = 0; stg < p.Stages; stg++ {
		/*@; BEGIN(CoA0=Nest) @*/
		for ic = 0; ic < p.Stages; ic++ {
			in.pa[stg][ic] = in.stepPrevious * p.a0[stg][ic]
		}

		stepStage := in.stepPrevious
		/*@; BEGIN(CoA1=Nest) @*/
		for ic = 0; ic < p.Stages; ic++ {
			stepStage *= in.stepRatio
			for id = 0; id < p.Stages; id++ {
				in.pa[stg][id] += p.cv[stg][ic] * stepStage * p.pv[ic][id]
			}
		}
	}
}
