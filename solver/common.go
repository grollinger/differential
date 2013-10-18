package solver

import "math"

type Function func(t float64, yT []float64, dy_out []float64)

type Config struct {
	// InitialStepSize, if > 0.0 specifies the step size
	// to be used in the first integration step
	// Else, the implementation should use a sensible default
	InitialStepSize float64

	// MinStepSize, if > 0.0 specifies the minimal size of a processing step
	// processing will abort, if this value could not be reached
	MinStepSize float64

	// MaxStepSize if > 0.0 specifies the maximum size of a processing step
	// processing will abort, if this value would be exceeded
	MaxStepSize float64

	//
	AbsoluteTolerance float64

	RelativeTolerance float64

	// MaxStepCount if > 0 specifies the maximum number number of steps the Integrator
	// will take before aborting processing if the target time has not been reached
	MaxStepCount uint

	// OneStepOnly, if set, causes the Integrator to stop processing
	// after the first integration step was performed
	OneStepOnly bool

	// Fcn contains the expression that should be evaluated for
	// the right hand side of the differential equation
	// yT'(t) = Fcn(t, yT(t))
	Fcn Function
}

type Statistics struct {
	// StepCount contains the number of steps the Integrator performed to calculate the
	StepCount uint
	// RejectedCount is the number of steps the Integrator rejected during processing
	RejectedCount uint
	// EvaluationCount is the number of times the right hand side expression
	// of the differential equation was evaluated during processing
	EvaluationCount uint

	// LastStepSize is the size of the last integration step performed
	LastStepSize float64
	// NextStepSize is the size of the next Step the integrator would take
	NextStepSize float64
	// CurrentTime is the value of t up to which the integration was performed
	CurrentTime float64
}

type Integrator interface {
	Info() IntegratorInfo
	Integrate(t, tEnd float64, yT []float64, config Config) (stat Statistics, err error)
}

type Problem interface {
	Initialize() []float64
	Fcn(t float64, yT []float64, dy_out []float64)
}

type IntegratorInfo struct {
	name          string
	stages, order uint
}

func (i *IntegratorInfo) Info() IntegratorInfo {
	return *i
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

func estimateStepSize(t float64, yT, fcnValue []float64, c *Config, order uint) float64 {
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
