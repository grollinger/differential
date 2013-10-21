package problems

import "testing"

func TestConstructor(t *testing.T) {
	b := NewBruss2D(2000).(*brusselator)
	if b.a != 3.4 {
		t.Error("Wrong Constant")
	}
}
