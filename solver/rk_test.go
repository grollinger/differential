package solver

import (
	"math"
	"math/rand"
	"testing"
)

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
	tMin, tMax float64
	sol        solution
	fcn        Function
	order      uint
	name       string
}

var integrationTests = []IntegrationTest{
	{-100, 100, order2, order2Deriv, 2, "x^2"},
	{-100, 100, order3, order3Deriv, 3, "x^3"},
	{-100, 100, order4, order4Deriv, 4, "x^4"},
}

func randomInInterval(low, high float64) float64 {
	return low + (rand.Float64() * (high - low))
}

func epsEqual(a, b, eps float64) bool {
	return math.Abs(a-b) < eps
}

func testIntegrators(t *testing.T) (integrators []Integrator) {
	integrators = make([]Integrator, NumberOfRKMethods+NumberOfPeerMethods)
	t.Logf("Integrators: %d", NumberOfRKMethods+NumberOfPeerMethods)
	var i uint = 0
	for j := 0; j < int(NumberOfRKMethods); j++ {
		rk, err := NewRK(RKMethod(j))
		if err != nil {
			t.Errorf("Couldn't create RK Method %d: %s", j, err.Error())
		} else {
			integrators[i] = rk
		}
		i++
	}
	for j := 0; j < int(NumberOfPeerMethods); j++ {
		p, err := NewPeer(PeerMethod(j))
		if err != nil {
			t.Errorf("Couldn't create Peer Method %d: %s", j, err.Error())
		} else {
			integrators[i] = p
		}
		t.Logf("Integrators #%d", i)
		i++
	}
	return
}

func TestRK(t *testing.T) {
	eps := 0.0001
	testsPerCase := 10

	methods := testIntegrators(t)
	for _, m := range methods {
		if m == nil {
			continue
		}

		info := m.Info()

		if testing.Verbose() {
			t.Logf("%s\tTest\tT0\tTE\tSteps\tReject\tEval\tLast h", info.name)
		}

		for _, v := range integrationTests {
			if v.order <= info.order {
				for i := 0; i < testsPerCase; i++ {
					t0 := randomInInterval(v.tMin, v.tMax)
					te := randomInInterval(t0, v.tMax)
					y := v.sol(t0)
					ye := v.sol(te)

					stat, err := m.Integrate(t0, te, y, Config{Fcn: v.fcn})

					if !epsEqual(stat.CurrentTime, te, eps) {
						t.Errorf("Tried to integrate up to %f but only reached %f", te, stat.CurrentTime)
					}
					if !epsEqual(y[0], ye[0], eps) {
						t.Errorf("Expected %f but result was %f", ye[0], y[0])
					}
					if err != nil {
						t.Errorf("Error: %s", err.Error())
					}
					if testing.Verbose() {
						t.Logf(" \t%s\t%.2f\t%.2f\t%d\t%d\t%d\t%.2f",
							v.name, t0, te, stat.StepCount, stat.RejectedCount, stat.EvaluationCount, stat.LastStepSize)
					}
				}
			} else {
				t.Logf("Skipped Test %s for RKMethod %s, order too high", v.name, info.name)
			}
		}
	}
}
