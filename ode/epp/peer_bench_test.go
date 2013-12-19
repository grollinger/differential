package epp

import (
	"fmt"
	. "github.com/rollingthunder/differential/ode"
	"github.com/rollingthunder/differential/problems"
	"testing"
)

func BenchmarkIntegrate(b *testing.B) {
	p, err := NewPeer(EPP4y3)

	if err != nil {
		b.Fatal(err)
	}

	brus := problems.NewBruss2D(200)
	y0 := brus.Initialize()

	cfg := Config{
		Fcn: brus.Fcn,
	}

	b.StartTimer()

	stat, err := p.Integrate(0.0, 1.0, y0, cfg)

	b.StopTimer()

	if testing.Verbose() {
		b.Log(fmt.Sprintln("Steps: %d, Evals: %d", stat.StepCount, stat.EvaluationCount))
	}
}
