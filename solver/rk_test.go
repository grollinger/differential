package solver

import "testing"

func TestRK2(t *testing.T) {
	_, err := NewRK(RK2)
	if err != nil {
		t.Errorf(err.Error())
	}
}
