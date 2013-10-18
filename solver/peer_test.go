package solver

import "testing"

func TestPeer(t *testing.T) {

	if !testing.Short() {
		t.Skipf("Skipping because we're NOT running in short test mode.")
	}

	peer, _ := NewPeer(EPP4)

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
