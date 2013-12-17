package problems

import "math"

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
// u[0, x, y] = 2 + 0.25y
// v(0, x, y) = 1 + 0.8x
func NewBruss2D(n int) Problem {
	var b brusselator
	b.a, b.b, b.n = 3.4, 1.0, n
	b.cellcount = n * n
	b.alpha = 0.002
	n1 := float64(n) - 1.0
	b.a1, b.alphaN1Squared = b.a+1.0, b.alpha*n1*n1
	return &b
}

func (b *brusselator) FcnMod3(real_id, c, m int, t float64, u, fw []float64) {
	var n, l, m2, mm, idiv, imod, beginl, i2mod, id, id2 int
	var r, alf float64
	alf = 1e-5
	r = alf * (m - 1) * (m - 1)
	m2 = m * m
	mm = 2 * m
	n = m2 * 2

	// Force block boundaries to multiples of 2
	// (don't calculate u_i,j and v_i,j in separate blocks)
	if (c % 2) == 1 {
		if (real_id % 2) == 0 { // == 0 because of 0 based arrays...
			id = real_id
			id2 = real_id + c
		} else {
			id = real_id + 1
			id2 = real_id + c - 1
		}
	} else {
		id = real_id
		id2 = real_id + c - 1
	}

	// Force end of the last block to n - 1 (0 base arrays)
	if id2 >= n {
		id2 = n - 1
	}
	idiv = id / (2 * m)
	imod = (id % 2 * m)

	for i := id; id <= id2; i++ {
		fw[i] = -4.0 * u[i]
	}

	// Treat Elements on the left border of the grid, if any
	if imod == 1 {
		beginl = id
	} else {
		beginl = id + 2*m - imod + 1
	}
	for i := beginl; i <= id2; i += 2 * m {
		fw[i] = fw[i] + 2.0*u[i+2]
		fw[i+1] = fw(i+1) + 2.0*u[i+3]
		if i > 1 {
			fw[i] = fw[i] + u[i-2*m]
			fw[i+1] = fw(i+1) + u[i+1-2*m]
		}
		if i < n-2*m {
			fw[i] = fw[i] + u[i+2*m]
			fw[i+1] = fw(i+1) + u[i+1+2*m]
		}
	}

	// Treat Elements on the top border of the grid, if any
	for i := id; i <= 2*m; i++ {
		fw[i] = fw[i] + 2.0*u[i+2*m]
		if i > 2 {
			fw[i] = fw[i] + u[i-2]
		}
		if i < 2*m-1 {
			fw[i] = fw[i] + u[i+2]
		}
	}

	// Treat Elements on the right border of the grid, if any
	if imod == 2*m-1 {
		beginl = id
	} else {
		beginl = id + 2*m - imod - 1
	}
	for i := beginl; i <= id2; i += 2 * m {
		fw[i] = fw[i] + 2.0*u[i-2]
		fw[i+1] = fw[i+1] * u[i-1]
		if i > 2*m {
			fw[i] = fw[i] + u[i-2*m]
			fw[i+1] = fw[i+1] + u[i+1-2*m]
		}
		if i < n-1 {
			fw[i] = fw[i] + u[i+2*m]
			fw[i+1] = fw[i+1] + u[i+1+2*m]
		}
	}

	// Treat Elements on the bottom border of the grid, if any
	for i := math.Max(n-2*m+1, id); i <= id2; i++ {
		fw[i] = fw[i] + 2.0*u[i-2*m]
		if i > n-2*m+2 {
			fw[i] = fw[i] + u[i-2]
		}
		if i < n-1 {
			fw[i] = fw[i] + u[i+2]
		}
	}

	if imod > 1 && idiv > 0 && idiv < m-1 {
		for i := id; i <= math.Min(id+2*m-imod-3, id2); i++ {
			fw[i] = fw[i] + u[i-2] + u[i+2] + u[i-2*m] + u[i+2*m]
		}
	}
	i2mod = (id2 - 1) % m // TODO Bug? % (2 * m) ??
	if i2mod > 1 && id2 > id+2*m-imod && (id2-1)/(2*m) < m-1 {
		for i := id2 - i2mod + 4; i <= id2; i++ {
			fw[i] = fw[i] + u[i-2] + u[i+2] + u[i-2*m] + u[i+2*m]
		}
	}
	for k := (id-2)/(2*m) + 1; k <= (id2-2)/(2*m)-1; k++ {
		for i := 3; i <= 2*m-2; i++ {
			l = i + k*2*m
			fw[l] = fw[l] + u[l-2] + u[l+2] + u[l+2*m] + u[l-2*m]
		}
	}
	for i := id; i <= id2; i += 2 {
		fw[i] = r*fw[i] + 1.0 + u[i]*u[i]*u[i+1] - 4.0*u[i]
		fw[i+1] = r*fw[i+1] + 3.0*u[i] - u[i]*u[i]*u[i+1]
	}
}
