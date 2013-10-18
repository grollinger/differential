package solver

import (
	"testing"
)

func TestRK(t *testing.T) {
	peer, _ := NewRK(DoPri5)

	testIntegrators(t, []Integrator{peer}, 1)
}
