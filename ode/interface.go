package ode

type Function func(t float64, yT []float64, dy_out []float64)
type BlockFunction func(startIdx, blockSize uint, t float64, yT []float64, dy_out []float64)

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

	// BlockSize if > 0 and <= n specifies how many entries of the system are evaluated in one go
	// must be <= n (the system size)
	// If 0, this will be forced to n
	// If FcnBlocked is unset, this will be forced to n
	BlockSize uint

	// Fcn or FcnBlocked contain the expression that should be evaluated for
	// the right hand side of the differential equation
	// yT'(t) = Fcn(t, yT(t))

	// Fcn
	Fcn        Function
	FcnBlocked BlockFunction
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
	Integrate(t, tEnd float64, yT []float64, config *Config) (stat Statistics, err error)
}

type IntegratorInfo struct {
	Name          string
	Stages, Order uint
}

func (i *IntegratorInfo) Info() IntegratorInfo {
	return *i
}
