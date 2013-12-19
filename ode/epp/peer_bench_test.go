package epp

import (
	. "github.com/rollingthunder/differential/ode"
	"github.com/rollingthunder/differential/problems"
	"testing"
)

func setupBenchmark() (p *peer, in integration, y0 []float64) {
	integrator, _ := NewPeer(EPP4y3)
	p = integrator.(*peer)

	brus := problems.NewBruss2D(200)
	y0 = brus.Initialize()

	cfg := Config{
		Fcn: brus.Fcn,
	}

	in = p.setupIntegration(0.0, 1.0, y0, cfg)

	return
}

func BenchmarkCoefficients(b *testing.B) {
	p, in, _ := setupBenchmark()

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		p.computeCoefficients(&in)
	}
}
