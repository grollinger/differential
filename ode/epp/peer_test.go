package epp

import (
	. "github.com/rollingthunder/differential/ode"
	. "github.com/rollingthunder/differential/ode/testing"
	"github.com/rollingthunder/differential/problems"
	"testing"
)

func TestPeer(t *testing.T) {
	peer, _ := NewPeer(EPP8_d)

	RunIntegratorTests(t, []Integrator{peer}, 1)
}

func TestAllPeer(t *testing.T) {
	if testing.Short() {
		t.Skipf("Skipping because we're running in short test mode.")
	}

	integrators := make([]Integrator, NumberOfPeerMethods)
	for j := 0; j < int(NumberOfPeerMethods); j++ {
		p, err := NewPeer(PeerMethod(j))
		if err != nil {
			t.Errorf("Couldn't create Peer Method %d: %s", j, err.Error())
		} else {
			integrators[j] = p
		}
	}

	RunIntegratorTests(t, integrators, 1)
}

func TestPeerBruss2D(t *testing.T) {
	peer, _ := NewPeer(EPP4)
	bruss := problems.NewBruss2D(10)
	instance := bruss.Initialize()

	config := Config{
		Fcn: bruss.Fcn,
	}

	peer.Integrate(0, 1, instance, config)
}

func TestPeerMBody4h(t *testing.T) {
	peer, _ := NewPeer(EPP8_d)
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
		t.Logf("MBody4H: %d steps, %d rejected, %d evaluations", stat.StepCount, stat.RejectedCount, stat.EvaluationCount)
		t.Logf("MBody: result[0..10] = %f", instance[:10])
	}
}
