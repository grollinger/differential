package problems

type cell struct {
	// Concentrations of the Reactants X and Y
	u, v float64
}

type brusselator struct {
	// Reaction Constants A and B
	a, b float64
	// Diffusion Constant alpha
	alpha float64
	// Grid Side length
	n int

	// Precalculated constants
	// alphaN1Squared = alpha * (n - 1)^2
	// a1 = A - 1
	alphaN1Squared, a1 float64

	concentration []cell
}

func v0(x, y float64) float64 {
	return 1 + 0.8*x
}

func u0(x, y float64) float64 {
	return 2 + 0.25*y
}

func (b *brusselator) Rate(index int) (du, dv float64) {
	top, right, bottom, left := index-b.n, index+1, index+b.n, index-1
	if top < 0 {
		top = bottom
	} else if bottom >= len(b.concentration) {
		bottom = top
	}

	if idxModN := index % b.n; idxModN == 0 {
		left = right
	} else if b.n-idxModN == 1 {
		right = left
	}

	conc := b.concentration

	du = b.b + conc[index].u*conc[index].u*conc[index].v - b.a1*conc[index].u + b.alphaN1Squared*(conc[top].u+conc[bottom].u+conc[left].u+conc[right].u-4.0*conc[index].u)

	dv = b.a*conc[index].u - conc[index].u*conc[index].u*conc[index].v + b.alphaN1Squared*(conc[top].v+conc[bottom].v+conc[left].v+conc[right].v-4.0*conc[index].v)

	return
}

// New initializes a new Brusselator 2D problem
// using the given grid side length n
// and standard parameters found in literature
// A = 3.4, B = 1.0
// alpha = 0.002
// u(0, x, y) = 2 + 0.25y
// v(0, x, y) = 1 + 0.8x
func NewBruss2D(n int) (b brusselator) {

	b.a, b.b, b.n = 3.4, 1.0, n
	b.alpha = 0.002
	n1 := float64(n) - 1.0
	b.a1, b.alphaN1Squared = b.a+1.0, b.alpha*n1*n1
	b.concentration = make([]cell, n*n)

	for i := range b.concentration {
		x, y := i%n, i/n
		xNorm, yNorm := float64(x)/n1, float64(y)/n1
		u, v := u0(xNorm, yNorm), v0(xNorm, yNorm)

		b.concentration[i] = cell{u, v}
	}

	return
}
