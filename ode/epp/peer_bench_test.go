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

	cfg.ValidateAndPrepare(uint(len(y0)))

	in = p.setupIntegration(0.0, 1.0, y0, &cfg)
	in.tCurrent, in.stepPrevious = p.startupIntegration(&in, 0.0)
	in.stepEstimate = in.stepPrevious
	return
}

func BenchmarkCoefficients(b *testing.B) {
	p, in, _ := setupBenchmark()

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		p.computeCoefficients(&in)
	}
}

func BenchmarkStages(b *testing.B) {
	p, in, _ := setupBenchmark()

	p.computeCoefficients(&in)

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		p.computeStages(&in)
	}
}

func BenchmarkEvaluations_NonBlocked(b *testing.B) {
	p, in, _ := setupBenchmark()

	p.computeCoefficients(&in)
	p.computeStages(&in)

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		p.computeEvaluations(&in)
	}
}

func BenchmarkErrorModel(b *testing.B) {
	p, in, _ := setupBenchmark()

	p.computeCoefficients(&in)
	p.computeStages(&in)
	p.computeEvaluations(&in)

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		p.computeErrorModel(&in)
	}
}
