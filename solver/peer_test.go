package solver

import "testing"

func TestPeer(t *testing.T) {
	peer, _ := NewPeer(EPP2)

	testIntegrators(t, []Integrator{peer}, 1)
}
