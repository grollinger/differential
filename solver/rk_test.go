package solver

import (
	"testing"
)

func TestRK(t *testing.T) {
	if !testing.Short() {
		t.Skipf("Skipping because we're NOT running in short test mode.")
	}

	peer, _ := NewRK(DoPri5)

	testIntegrators(t, []Integrator{peer}, 1)
}

func TestAllRK(t *testing.T) {
	if testing.Short() {
		t.Skipf("Skipping because we're running in short test mode.")
	}

	integrators := make([]Integrator, NumberOfRKMethods+NumberOfPeerMethods)
	for j := 0; j < int(NumberOfRKMethods); j++ {
		rk, err := NewRK(RKMethod(j))
		if err != nil {
			t.Errorf("Couldn't create RK Method %d: %s", j, err.Error())
		} else {
			integrators[j] = rk
		}
	}

	testIntegrators(t, integrators, iterationsPerTest)
}
