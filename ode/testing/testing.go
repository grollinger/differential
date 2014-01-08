package testing

import (
	. "github.com/rollingthunder/differential/ode"
	"github.com/rollingthunder/differential/util"
	"math"
	"testing"
)

var iterationsPerTest = 10

type solution func(float64) []float64

// x^2 + C
func order2(t float64) []float64 {
	return []float64{t * t}
}
func order2Deriv(t float64, y []float64, dy []float64) {
	dy[0] = 2 * t
}

// 10*x^3+ PI * x^2 +142 * x + 10
func order3(t float64) []float64 {
	return []float64{10*math.Pow(t, 3) + math.Pi*math.Pow(t, 2) + 142*t + 10}
}
func order3Deriv(t float64, y []float64, dy []float64) {
	dy[0] = 30*math.Pow(t, 2) + 2*math.Pi*t + 142
}

// 83.23454 * x^4 + 42.4543*x^3+ E * x^2
func order4(t float64) []float64 {
	return []float64{83.23454*math.Pow(t, 4) + 42.4543*math.Pow(t, 3) + math.E*math.Pow(t, 2)}
}
func order4Deriv(t float64, y []float64, dy []float64) {
	dy[0] = 4*83.23454*math.Pow(t, 3) + 3*42.4543*math.Pow(t, 2) + 2*math.E*t
}

type IntegrationTest struct {
	TMin, TMax float64
	Sol        solution
	Fcn        Function
	Order      uint
	Name       string
}

var IntegrationTests = []IntegrationTest{
	{-10, 10, order2, order2Deriv, 2, "x^2"},
	{-10, 10, order3, order3Deriv, 3, "x^3"},
	{-10, 10, order4, order4Deriv, 4, "x^4"},
}

func RunIntegratorTests(t *testing.T, methods []Integrator, iterations int) {
	var eps float64 = 0.0001

	for _, m := range methods {
		if m == nil {
			continue
		}

		info := m.Info()

		if testing.Verbose() {
			t.Logf("%s\tTest\tT0\tTE\tSteps\tReject\tEval\tLast h", info.Name)
		}

		for _, v := range IntegrationTests {
			if v.Order <= info.Order {
				for i := 0; i < iterations; i++ {
					t0 := util.RandomInInterval(v.TMin, v.TMax)
					te := util.RandomInInterval(t0, v.TMax)
					y := v.Sol(t0)
					ye := v.Sol(te)

					stat, err := m.Integrate(t0, te, y, &Config{Fcn: v.Fcn})

					if !util.EpsEqual(stat.CurrentTime, te, eps) {
						t.Errorf("Tried to integrate up to %f but only reached %f", te, stat.CurrentTime)
					}
					if !util.EpsEqual(y[0], ye[0], eps) {
						t.Errorf("Expected %f but result was %f", ye[0], y[0])
					}
					if err != nil {
						t.Errorf("Error: %s", err.Error())
					}
					if testing.Verbose() {
						t.Logf(" \t%s\t%.2f\t%.2f\t%d\t%d\t%d\t%.2f",
							v.Name, t0, te, stat.StepCount, stat.RejectedCount, stat.EvaluationCount, stat.LastStepSize)
					}
				}
			} else {
				t.Logf("Skipped Test %s for RKMethod %s, order too high", v.Name, info.Name)
			}
		}
	}
}
