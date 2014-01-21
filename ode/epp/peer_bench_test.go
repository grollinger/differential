package epp

import (
	"fmt"
	. "github.com/rollingthunder/differential/ode"
	"github.com/rollingthunder/differential/problems"
	"github.com/rollingthunder/differential/util"
	"strconv"
	"testing"
)

func setupBruss() (p *peer, in integration, y0 []float64) {
	return setupBenchmark(problems.NewBruss2D(200))
}

func setupBenchmark(prob problems.TiledProblem) (p *peer, in integration, y0 []float64) {
	integrator, _ := NewPeer(EPP4y3)
	p = integrator.(*peer)

	y0 = prob.Initialize()

	cfg := Config{
		Fcn: prob.Fcn,
	}

	cfg.ValidateAndPrepare(uint(len(y0)))

	in = p.setupIntegration(0.0, 1.0, y0, &cfg)
	in.tCurrent, in.stepPrevious = p.startupIntegration(&in, 0.0)
	in.stepEstimate = in.stepPrevious
	return
}

func BenchmarkCoefficients(b *testing.B) {
	p, in, _ := setupBruss()

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		p.computeCoefficients(&in)
	}
}

var problemVariants = []struct {
	Name        string
	Constructor func(uint) problems.TiledProblem
}{
	{"Bruss2D", problems.NewBruss2D},
}

var sizeVariants = []uint{
	5,
	10,
	15,
	30,
}

type namedImplementation struct {
	Name string
	Impl computationStep
}

func benchmarkComputationStep(stepName string, prepareIntegration computationStep, implementations []namedImplementation) {
	const TIME string = "Time"
	const NORMALIZED_TIME string = "Time/n"

	var resultTables = make([]util.Table, 0, len(problemVariants))

	for i_problem := range problemVariants {
		prob := problemVariants[i_problem]
		var result = util.Table{
			Title:      fmt.Sprintf("%v - %v", stepName, prob.Name),
			RowHeaders: make([]string, 0, len(implementations)),
			ColHeaders: make([]string, 0, len(sizeVariants)),
			Data: map[string][][]float64{
				TIME:            util.MakeRectangular(uint(len(implementations)), uint(len(sizeVariants))),
				NORMALIZED_TIME: util.MakeRectangular(uint(len(implementations)), uint(len(sizeVariants))),
			},
		}
		for j_size := range sizeVariants {
			s := sizeVariants[j_size]
			pConfig := prob.Constructor(s)

			for k_impl := range implementations {
				variant := implementations[k_impl]
				p, in, _ := setupBenchmark(pConfig)

				prepareIntegration(p, &in)

				if k_impl == 0 {
					result.ColHeaders = append(result.ColHeaders, strconv.Itoa(int(in.n)))
				}

				if j_size == 0 {
					result.RowHeaders = append(result.RowHeaders, variant.Name)
				}

				benchResult := testing.Benchmark(
					func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							variant.Impl(p, &in)
						}
					})
				time := benchResult.NsPerOp()

				result.Data[TIME][k_impl][j_size] = float64(time)
				result.Data[NORMALIZED_TIME][k_impl][j_size] = float64(time) / float64(in.n)
			}
		}
		resultTables = append(resultTables, result)
	}
	util.WriteTablesFile(resultTables, fmt.Sprintf("%s.html", stepName))
}

func TestBenchmarkStages(t *testing.T) {
	var prepare computationStep = func(p *peer, in *integration) {
		p.computeStages(in)
	}

	var stagesVariants = []namedImplementation{
		{"Vanilla", (*peer).computeStages},
		{"FuseAB", (*peer).computeStages_FuseAB},
		{"FuseAB_ExchangeIJ", (*peer).computeStages_FuseAB_ExchangeIJ},
	}

	benchmarkComputationStep(
		"Stages",
		prepare,
		stagesVariants,
	)
}

func BenchmarkEvaluations_NonBlocked(b *testing.B) {
	p, in, _ := setupBruss()

	p.computeCoefficients(&in)
	p.computeStages(&in)

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		p.computeEvaluations(&in)
	}
}

func BenchmarkErrorModel(b *testing.B) {
	p, in, _ := setupBruss()

	p.computeCoefficients(&in)
	p.computeStages(&in)
	p.computeEvaluations(&in)

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		p.computeErrorModel(&in)
	}
}
