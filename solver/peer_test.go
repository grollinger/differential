package solver

import (
	"github.com/rollingthunder/differential/problems"
	"testing"
)

func TestPeer(t *testing.T) {
	peer, _ := NewPeer(EPP2)

	testIntegrators(t, []Integrator{peer}, 1)
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

	testIntegrators(t, integrators, iterationsPerTest)
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
