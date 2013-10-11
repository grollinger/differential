package solver

type peer struct {
	parameters
	method        PeerMethod
	stages, order uint
	icmin, icmax  int

	sigx, ema, emaoh float64
	c, e             []float64
	b, a0, cv, pv    [][]float64
}
