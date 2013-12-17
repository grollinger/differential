package problems

import "github.com/rollingthunder/differential/util"

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
// using the given grid side length N
// and standard parameters found in literature
// A = 3.4, B = 1.0
// alpha = 0.002
// u[0, x, y] = 2 + 0.25y
// v(0, x, y) = 1 + 0.8x
func NewBruss2D(N int) Problem {
	var b brusselator
	b.a, b.b, b.n = 3.4, 1.0, N
	b.cellcount = N * N
	b.alpha = 0.002
	n1 := float64(N) - 1.0
	b.a1, b.alphaN1Squared = b.a+1.0, b.alpha*n1*n1
	return &b
}

func (b *brusselator) FcnBlock(startIdx, blockSize int, t float64, yT []float64, dy_out []float64) {
	var l, idiv, imod, beginl, i2mod, id, id2 int

	systemSize := b.cellcount * 2

	// Force block boundaries to multiples of 2
	// (don't calculate u_i,j and v_i,j in separate blocks)
	if (blockSize % 2) == 1 {
		if (startIdx % 2) == 0 { // == 0 because of 0 based arrays...
			id = startIdx
			id2 = startIdx + blockSize
		} else {
			id = startIdx + 1
			id2 = startIdx + blockSize - 1
		}
	} else {
		id = startIdx
		id2 = startIdx + blockSize - 1
	}

	// Force end of the last block to
	if id2 >= systemSize {
		id2 = systemSize - 2
	}
	idiv = id / (2 * b.n)
	imod = id % (2 * b.n)

	for i := id; i <= id2; i++ {
		dy_out[i] = -4.0 * yT[i]
	}

	// Treat Elements on the left border of the grid, if any
	if imod == 0 {
		beginl = id
	} else {
		beginl = id + 2*b.n - imod
	}
	for i := beginl; i <= id2; i += 2 * b.n {
		dy_out[i] = dy_out[i] + 2.0*yT[i+2]
		dy_out[i+1] = dy_out[i+1] + 2.0*yT[i+3]
		if i > 1 && i < systemSize-2*b.n {
			dy_out[i] = dy_out[i] + yT[i-2*b.n]
			dy_out[i+1] = dy_out[i+1] + yT[i+1-2*b.n]
			dy_out[i] = dy_out[i] + yT[i+2*b.n]
			dy_out[i+1] = dy_out[i+1] + yT[i+1+2*b.n]
		}
	}

	// Treat Elements on the top border of the grid, if any
	for i := id; i <= 2*b.n-1; i++ {
		dy_out[i] = dy_out[i] + 2.0*yT[i+2*b.n]
		// only for Elements NOT on the right or left borders
		// those cases are handled elsewhere
		if i > 1 && i < 2*b.n-2 {
			dy_out[i] = dy_out[i] + yT[i-2]
			dy_out[i] = dy_out[i] + yT[i+2]
		}
	}

	// Treat Elements on the right border of the grid, if any
	endl := -1
	if imod == 2*b.n-2 {
		endl = id
	} else {
		endl = id + 2*b.n - imod - 2
	}
	for i := endl; i <= id2; i += 2 * b.n {
		dy_out[i] = dy_out[i] + 2.0*yT[i-2]
		dy_out[i+1] = dy_out[i+1] + 2.0*yT[i-1]
		if i >= 2*b.n && i < systemSize-2*b.n {
			dy_out[i] = dy_out[i] + yT[i-2*b.n]
			dy_out[i+1] = dy_out[i+1] + yT[i+1-2*b.n]
			dy_out[i] = dy_out[i] + yT[i+2*b.n]
			dy_out[i+1] = dy_out[i+1] + yT[i+1+2*b.n]
		}
	}

	// Treat Elements on the bottom border of the grid, if any
	for i := util.Max(systemSize-2*b.n, id); i <= id2; i++ {
		dy_out[i] = dy_out[i] + 2.0*yT[i-2*b.n]
		// only for Elements NOT on the right or left borders
		// those cases are handled elsewhere
		if i > systemSize-2*b.n+1 && i < systemSize-2 {
			dy_out[i] = dy_out[i] + yT[i-2]
			dy_out[i] = dy_out[i] + yT[i+2]
		}
	}
	// Elements not in the top or bottom row, but in the first
	// and potentially incomplete row of the block
	if imod > 0 && idiv > 0 && idiv < b.n-1 {
		for i := id; i <= util.Min(id+2*b.n-imod-2, id2); i++ {
			dy_out[i] = dy_out[i] + yT[i-2] + yT[i+2] + yT[i-2*b.n] + yT[i+2*b.n]
		}
	}
	i2mod = (id2 - 1) % (2 * b.n)
	if i2mod > 1 && id2 >= id+2*b.n-imod && (id2)/(2*b.n) < b.n-1 {
		for i := id2 - i2mod + 1; i <= id2; i++ {
			dy_out[i] = dy_out[i] + yT[i-2] + yT[i+2] + yT[i-2*b.n] + yT[i+2*b.n]
		}
	}
	for k := util.Max((id-1)/(2*b.n)+1, 1); k <= util.Min((id2+1)/(2*b.n)-1, b.n-2); k++ {
		for i := 2; i < 2*b.n-2; i++ {
			l = i + k*2*b.n
			dy_out[l] = dy_out[l] + yT[l-2] + yT[l+2] + yT[l+2*b.n] + yT[l-2*b.n]
		}
	}
	for i := id; i <= id2; i += 2 {
		dy_out[i] = b.alphaN1Squared*dy_out[i] + b.b + yT[i]*yT[i]*yT[i+1] - b.a1*yT[i]
		dy_out[i+1] = b.alphaN1Squared*dy_out[i+1] + b.a*yT[i] - yT[i]*yT[i]*yT[i+1]
	}
}
