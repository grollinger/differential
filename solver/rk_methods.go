package solver

import "errors"

const (
	RK2               = RKMethod(iota) // RK2(3)
	RKFB4                              // RKFB4(5)
	DoPri5                             // DoPri5(4)
	NumberOfRKMethods = uint(iota)
)

func NewRK(m RKMethod) (i Integrator, err error) {
	var r rk
	switch m {
	case RK2:
		r.stages, r.order = 3, 3
		r.name = "RK2"
		makeCoeffs(&r)
		setCoeffsRK2(&r)
	case RKFB4:
		r.stages, r.order = 6, 4
		r.name = "RKFB4"
		makeCoeffs(&r)
		setCoeffsRKFB4(&r)
	case DoPri5:
		r.stages, r.order = 7, 5
		r.name = "DoPri5"
		r.firstStageAsLast = true
		makeCoeffs(&r)
		setCoeffsDoPri5(&r)

	default:
		err = errors.New("unknown rk method")
	}

	i = &r
	return
}

func makeCoeffs(r *rk) {
	r.b, r.c, r.e = make([]float64, r.stages), make([]float64, r.stages), make([]float64, r.stages)
	r.a = makeSquare(r.stages)
}
func setCoeffsRK2(r *rk) {
	r.a[1][0] = 1.0
	r.a[2][0] = 1.0 / 4.0
	r.a[2][1] = 1.0 / 4.0

	r.c[0] = 0.0
	r.c[1] = 1.0
	r.c[2] = 0.5
	r.b[0] = 1.0 / 6.0
	r.b[1] = 1.0 / 6.0
	r.b[2] = 2.0 / 3.0
	r.e[0] = -1.0 / 3.0
	r.e[1] = -1.0 / 3.0
	r.e[2] = 2.0 / 3.0
}
func setCoeffsRKFB4(r *rk) {
	r.a[1][0] = 1.0 / 4
	r.a[2][0] = 3.0 / 32.0
	r.a[2][1] = 9.0 / 32.0
	r.a[3][0] = 1932.0 / 2197.0
	r.a[3][1] = -7200.0 / 2197.0
	r.a[3][2] = 7296.0 / 2197.0
	r.a[4][0] = 439.0 / 216.0
	r.a[4][1] = -8.0
	r.a[4][2] = 3680.0 / 513.0
	r.a[4][3] = -845.0 / 4104.0
	r.a[5][0] = -8.0 / 27.0
	r.a[5][1] = 2.0
	r.a[5][2] = -3544.0 / 2565.0
	r.a[5][3] = 1859.0 / 4104.0
	r.a[5][4] = -11.0 / 40.0

	r.c[0] = 0.0
	r.c[1] = 1.0 / 4.0
	r.c[2] = 3.0 / 8
	r.c[3] = 12.0 / 13.0
	r.c[4] = 1.0
	r.c[5] = 1.0 / 2.0

	r.b[0] = 25.0 / 216.0
	r.b[1] = 0.0
	r.b[2] = 1408.0 / 2565.0
	r.b[3] = 2197.0 / 4104.0
	r.b[4] = -1.0 / 5.0
	r.b[5] = 0.0

	r.e[0] = 16.0 / 135.0
	r.e[1] = 0.0
	r.e[2] = 6656.0 / 12825.0
	r.e[3] = 28561.0 / 56430.0
	r.e[4] = -9.0 / 50.0
	r.e[5] = 2.0 / 55.0

	// for difference of solutions
	var i uint
	for i = 0; i < r.stages; i++ {
		r.e[i] = r.b[i] - r.e[i]
	}
}
func setCoeffsDoPri5(r *rk) {
	r.a[1][0] = 0.2
	r.a[2][0] = 3.0 / 40.0
	r.a[2][1] = 9.0 / 40.0
	r.a[3][0] = 44.0 / 45.0
	r.a[3][1] = -56.0 / 15.0
	r.a[3][2] = 32.0 / 9.0
	r.a[4][0] = 19372.0 / 6561.0
	r.a[4][1] = -25360.0 / 2187.0
	r.a[4][2] = 64448.0 / 6561.0
	r.a[4][3] = -212.0 / 729.0
	r.a[5][0] = 9017.0 / 3168.0
	r.a[5][1] = -355.0 / 33.0
	r.a[5][2] = 46732.0 / 5247.0
	r.a[5][3] = 49.0 / 176.0
	r.a[5][4] = -5103.0 / 18656.0
	r.a[6][0] = 35.0 / 384.0
	r.a[6][1] = 0.0
	r.a[6][2] = 500.0 / 1113.0
	r.a[6][3] = 125.0 / 192.0
	r.a[6][4] = -2187.0 / 6784.0
	r.a[6][5] = 11.0 / 84.0

	r.b[0] = 35.0 / 384.0
	r.b[1] = 0.0
	r.b[2] = 500.0 / 1113.0
	r.b[3] = 125.0 / 192.0
	r.b[4] = -2187.0 / 6784.0
	r.b[5] = 11.0 / 84.0
	r.b[6] = 0.0

	r.c[0] = 0.0
	r.c[1] = 0.2
	r.c[2] = 0.3
	r.c[3] = 0.8
	r.c[4] = 8.0 / 9.0
	r.c[5] = 1.0
	r.c[6] = 1.0

	r.e[0] = 71.0 / 57600.0
	r.e[1] = 0.0
	r.e[2] = -71.0 / 16695.0
	r.e[3] = 71.0 / 1920.0
	r.e[4] = -17253.0 / 339200.0
	r.e[5] = 22.0 / 525.0
	r.e[6] = -1.0 / 40.0
}
