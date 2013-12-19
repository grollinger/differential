package rk

import (
	. "github.com/rollingthunder/differential/ode"
	. "github.com/rollingthunder/differential/ode/testing"
	"github.com/rollingthunder/differential/problems"
	"testing"
)

func TestAllRK(t *testing.T) {
	if testing.Short() {
		t.Skipf("Skipping because we're running in short test mode.")
	}

	integrators := make([]Integrator, NumberOfRKMethods)
	for j := 0; j < int(NumberOfRKMethods); j++ {
		rk, err := NewRK(RKMethod(j))
		if err != nil {
			t.Errorf("Couldn't create RK Method %d: %s", j, err.Error())
		} else {
			integrators[j] = rk
		}
	}

	RunIntegratorTests(t, integrators, 1)
}

func TestRKMBody4h(t *testing.T) {
	peer, _ := NewRK(DoPri5)
	mbody := problems.NewMBody(4)
	instance := mbody.Initialize()

	config := Config{
		Fcn:               mbody.Fcn,
		AbsoluteTolerance: 1.e-5,
		RelativeTolerance: 1.e-5,
	}
	t0, te := 0.0, 0.1

	stat, err := peer.Integrate(t0, te, instance, config)

	if err != nil {
		t.Fatalf("Integration failed - %s", err.Error())
	}

	if testing.Verbose() {
		t.Logf("MBody: %d steps, %d rejected, %d evaluations", stat.StepCount, stat.RejectedCount, stat.EvaluationCount)
		t.Logf("MBody: result[0..10] = %f", instance[:10])
	}
}
