package problems

type brusselator struct {
	// Reaction Constants A and B
	a, b float64
	// Diffusion Constant alpha
	alpha float64
	// Grid Side length
	n int
	// Grid Cell Count for input validation
	cellcount int

	// Precalculated constants
	// alphaN1Squared = alpha * (n - 1)^2
	// a1 = A - 1
	alphaN1Squared, a1 float64
}

func v0(x, y float64) float64 {
	return 1 + 0.8*x
}

func u0(x, y float64) float64 {
	return 2 + 0.25*y
}

func (b *brusselator) Initialize() (grid []float64) {
	grid = make([]float64, b.cellcount*2)

	n1 := float64(b.n) - 1.0

	for i := range grid {
		if i%2 == 1 {
			continue
		}

		x, y := i%b.n, i/b.n
		xNorm, yNorm := float64(x)/n1, float64(y)/n1
		u, v := u0(xNorm, yNorm), v0(xNorm, yNorm)

		grid[i] = u
		grid[i+1] = v
	}
	return
}

func (b *brusselator) Fcn(t float64, yT []float64, dy_out []float64) {
	// TODO validate input

	for index := 0; index < b.cellcount; index++ {
		top, right, bottom, left := index-b.n, index+1, index+b.n, index-1
		if top < 0 {
			top = bottom
		} else if bottom >= b.cellcount {
			bottom = top
		}

		if idxModN := index % b.n; idxModN == 0 {
			left = right
		} else if b.n-idxModN == 1 {
			right = left
		}
		top <<= 1
		bottom <<= 1
		left <<= 1
		right <<= 1
		here := index << 1

		// du = B + u^2*v - (A + 1)*u + (alpha * (n-1)^2) (u_top+u_bottom+u_left+u_right-4 * u)
		dy_out[here] = b.b + yT[here]*yT[here]*yT[here+1] - b.a1*yT[here] + b.alphaN1Squared*(yT[top]+yT[bottom]+yT[left]+yT[right]-4.0*yT[here])

		// dv = A * u - u^2*v + (alpha * (n-1)^2) * (v_top + v_bottom + v_left + v_right - 4 * v)
		dy_out[here+1] = b.a*yT[here] - yT[here]*yT[here]*yT[here+1] + b.alphaN1Squared*(yT[top+1]+yT[bottom+1]+yT[left+1]+yT[right+1]-4.0*yT[here+1])

	}
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
	b.cellcount = n * n
	b.alpha = 0.002
	n1 := float64(n) - 1.0
	b.a1, b.alphaN1Squared = b.a+1.0, b.alpha*n1*n1
	return
}
