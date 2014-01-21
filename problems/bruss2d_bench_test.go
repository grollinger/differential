package problems

import "testing"

var benchmarkSystemSize uint = 30
var benchmarkT float64 = 0.0

func setupBenchmark() (b *brusselator, instance, out []float64) {
	b = NewBruss2D(benchmarkSystemSize).(*brusselator)
	instance = b.Initialize()
	out = make([]float64, len(instance))
	return
}

func BenchmarkCorrect(b *testing.B) {
	brus, inst, out := setupBenchmark()

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		brus.FcnCorrect(0, inst, out)
	}
}

func benchmarkBlocked(blockSize uint, b *testing.B) {
	brus, inst, out := setupBenchmark()

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		var block uint
		for block = 0; block < uint(len(inst)); block += blockSize {
			brus.FcnBlock(block, blockSize, benchmarkT, inst, out)
		}
	}
}

func BenchmarkBlocked1(b *testing.B)   { benchmarkBlocked(1, b) }
func BenchmarkBlocked2(b *testing.B)   { benchmarkBlocked(2, b) }
func BenchmarkBlocked3(b *testing.B)   { benchmarkBlocked(3, b) }
func BenchmarkBlocked4(b *testing.B)   { benchmarkBlocked(4, b) }
func BenchmarkBlocked5(b *testing.B)   { benchmarkBlocked(5, b) }
func BenchmarkBlocked10(b *testing.B)  { benchmarkBlocked(10, b) }
func BenchmarkBlocked20(b *testing.B)  { benchmarkBlocked(20, b) }
func BenchmarkBlocked40(b *testing.B)  { benchmarkBlocked(40, b) }
func BenchmarkBlocked50(b *testing.B)  { benchmarkBlocked(50, b) }
func BenchmarkBlocked60(b *testing.B)  { benchmarkBlocked(60, b) }
func BenchmarkBlocked100(b *testing.B) { benchmarkBlocked(100, b) }
