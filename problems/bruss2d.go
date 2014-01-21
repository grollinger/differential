package problems

import "github.com/rollingthunder/differential/util"
import "fmt"

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

func (b *brusselator) Description() string {
	return fmt.Sprintf("Brusselator2D with %d x %d cells", b.n, b.n)
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

func (b *brusselator) FcnCorrect(t float64, yT []float64, dy_out []float64) {
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
func NewBruss2D(N uint) TiledProblem {
	if N <= 0 {
		return nil
	}

	var b brusselator
	b.a, b.b, b.n = 3.4, 1.0, int(N)
	b.cellcount = int(N * N)
	b.alpha = 0.002
	n1 := float64(N) - 1.0
	b.a1, b.alphaN1Squared = b.a+1.0, b.alpha*n1*n1
	return &b
}

func (b *brusselator) Fcn(t float64, yT []float64, dy_out []float64) {
	b.FcnBlock(0, uint(len(yT)), t, yT, dy_out)
}

func (b *brusselator) FcnBlock(startIdx, blockSize uint, t float64, yT []float64, dy_out []float64) {
	var lo, hi int
	var loX, loY, hiX int

	systemSize := b.cellcount * 2
	line := 2 * b.n

	// Force block boundaries to multiples of 2
	// (don't calculate u_i,j and v_i,j in separate blocks)
	if (blockSize % 2) == 1 {
		if (startIdx % 2) == 0 { // == 0 because of 0 based arrays...
			lo = int(startIdx)
			hi = int(startIdx + blockSize)
		} else {
			lo = int(startIdx + 1)
			hi = int(startIdx + blockSize - 1)
		}
	} else {
		lo = int(startIdx)
		hi = int(startIdx + blockSize - 1)
	}

	// Force end of the last block to
	if hi >= systemSize {
		hi = systemSize - 1
	}
	loY = lo / line
	loX = lo % line
	hiX = (hi - 1) % line

	for i := lo; i <= hi; i++ {
		dy_out[i] = -4.0 * yT[i]
	}

	// Treat Elements on the left border of the grid, if any
	var lineBeginning int
	if loX == 0 {
		lineBeginning = lo
	} else {
		lineBeginning = lo + line - loX
	}
	for i := lineBeginning; i <= hi; i += line {
		dy_out[i] += 2.0 * yT[i+2]
		dy_out[i+1] += 2.0 * yT[i+3]
		if i > 1 && i < systemSize-line {
			dy_out[i] += yT[i-line]
			dy_out[i+1] += yT[i+1-line]
			dy_out[i] += yT[i+line]
			dy_out[i+1] += yT[i+1+line]
		}
	}

	// Treat Elements on the top border of the grid, if any
	for i := lo; i <= line-1; i++ {
		dy_out[i] = dy_out[i] + 2.0*yT[i+line]
		// only for Elements NOT on the right or left borders
		// those cases are handled elsewhere
		if i > 1 && i < line-2 {
			dy_out[i] += yT[i-2]
			dy_out[i] += yT[i+2]
		}
	}

	// Treat Elements on the right border of the grid, if any
	var lineEnd int
	if loX == line-2 {
		lineEnd = lo
	} else {
		lineEnd = lo + line - loX - 2
	}
	for i := lineEnd; i <= hi; i += line {
		dy_out[i] += 2.0 * yT[i-2]
		dy_out[i+1] += 2.0 * yT[i-1]
		if i >= line && i < systemSize-line {
			dy_out[i] += yT[i-line]
			dy_out[i+1] += yT[i+1-line]
			dy_out[i] += yT[i+line]
			dy_out[i+1] += yT[i+1+line]
		}
	}

	// Treat Elements on the bottom border of the grid, if any
	for i := util.Max(systemSize-line, lo); i <= hi; i++ {
		dy_out[i] += 2.0 * yT[i-line]
		// only for Elements NOT on the right or left borders
		// those cases are handled elsewhere
		if i > systemSize-line+1 && i < systemSize-2 {
			dy_out[i] += yT[i-2]
			dy_out[i] += yT[i+2]
		}
	}
	// Elements not in the top or bottom row, but in the first
	// and potentially incomplete row of the block
	if loX > 0 && loY > 0 && loY < b.n-1 {
		for i := lo; i <= util.Min(lineEnd-1, hi); i++ {
			dy_out[i] += yT[i-2] + yT[i+2] + yT[i-line] + yT[i+line]
		}
	}

	if hiX > 1 && hi >= lo+line-loX && (hi)/(line) < b.n-1 {
		for i := hi - hiX + 1; i <= hi; i++ {
			dy_out[i] = dy_out[i] + yT[i-2] + yT[i+2] + yT[i-line] + yT[i+line]
		}
	}
	for k := util.Max((lo-1)/line+1, 1); k <= util.Min((hi+1)/line-1, b.n-2); k++ {
		for i := 2; i < line-2; i++ {
			l := i + k*line
			dy_out[l] = dy_out[l] + yT[l-2] + yT[l+2] + yT[l+line] + yT[l-line]
		}
	}
	for i := lo; i <= hi; i += 2 {
		dy_out[i] = b.alphaN1Squared*dy_out[i] + b.b + yT[i]*yT[i]*yT[i+1] - b.a1*yT[i]
		dy_out[i+1] = b.alphaN1Squared*dy_out[i+1] + b.a*yT[i] - yT[i]*yT[i]*yT[i+1]
	}
}
