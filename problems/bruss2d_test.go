package problems

import "testing"
import "github.com/rollingthunder/differential/util"

func TestConstructor(t *testing.T) {
	b := NewBruss2D(2000).(*brusselator)
	if b.a != 3.4 {
		t.Error("Wrong Constant")
	}
}

func TestBlock_even(t *testing.T) {
	blockSize := 100
	b := NewBruss2D(50).(*brusselator)
	var yT = b.Initialize()
	var block, nonBlock = make([]float64, len(yT)), make([]float64, len(yT))
	b.Fcn(0, yT, nonBlock)
	for i := 0; i < len(yT); i += blockSize {
		b.FcnBlock(i, blockSize, 0.0, yT, block)
	}
	t.Log(nonBlock[0:5])
	t.Log(block[0:5])

	if !util.ArrayEpsEquals(block, nonBlock, 1e-3) {
		t.Errorf("Blocking and non blocking variants have different results")
	}
}

func TestBlock_odd(t *testing.T) {
	blockSize := 101
	b := NewBruss2D(50).(*brusselator)
	var yT = b.Initialize()
	var block, nonBlock = make([]float64, len(yT)), make([]float64, len(yT))
	b.Fcn(0, yT, nonBlock)
	for i := 0; i < len(yT); i += blockSize {
		b.FcnBlock(i, blockSize, 0.0, yT, block)
	}
	t.Log(nonBlock[0:5])
	t.Log(block[0:5])

	if !util.ArrayEpsEquals(block, nonBlock, 1e-3) {
		t.Errorf("Blocking and non blocking variants have different results")
	}
}
