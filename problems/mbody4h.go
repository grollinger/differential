package problems

import "math"

const mbody = 4
const meps = 1e-4

type mbody4h struct {
	mass []float64
}

func NewMBody4h(n uint) (p Problem) {
	var m mbody4h

	m.mass = make([]float64, mbody)

	rf1 := 4 * math.Pi / 8

	for i := range m.mass {
		ip := float64(i + 1)
		m.mass[i] = (0.3 + 0.1*(math.Cos(ip*rf1)+1.0)) / mbody
	}

	return &m
}

func (m *mbody4h) Initialize() (y0 []float64) {
	y0 = make([]float64, mbody*6)

	rf1 := 4 * math.Pi / 8
	rf2 := 2 * math.Pi / mbody

	for i := range m.mass {
		i1 := float64(i + 1)
		m.mass[i] = (0.3 + 0.1*(math.Cos(i1*rf1)+1.0)) / mbody
		rad := 1.7 + math.Cos(i1*0.75)
		v := 0.22 * math.Sqrt(rad)
		ci := math.Cos(i1 * rf2)
		si := math.Sin(i1 * rf2)

		ip := 6 * i
		y0[ip] = rad * ci
		y0[ip+1] = rad * si
		y0[ip+2] = 0.4 * si
		y0[ip+3] = -v * si
		y0[ip+4] = v * ci
		y0[ip+5] = 0
	}
	return
}

func (m *mbody4h) Fcn(t float64, yT []float64, dy_out []float64) {

	for i := 0; i < mbody; i++ {
		ip := 6 * i
		dy_out[ip] = yT[ip+3]
		dy_out[ip+1] = yT[ip+4]
		dy_out[ip+2] = yT[ip+5]
		var f1, f2, f3 float64 = 0, 0, 0

		for j := 0; j < mbody; j++ {
			if i != j {
				jp := 6 * j
				dist := meps + math.Pow(yT[ip]-yT[jp], 2) + math.Pow(yT[ip+1]-yT[jp+1], 2) + math.Pow(yT[ip+2]-yT[jp+2], 2)
				dist = m.mass[j] / (dist * math.Sqrt(dist))
				f1 = f1 + (yT[jp]-yT[ip])*dist
				f2 = f2 + (yT[jp+1]-yT[ip+1])*dist
				f3 = f3 + (yT[jp+2]-yT[ip+2])*dist
			}
		}

		dy_out[ip+3] = f1
		dy_out[ip+4] = f2
		dy_out[ip+5] = f3
	}
}

func (m *mbody4h) invariant(n int, yT []float64) float64 {

	en := 0.0

	for i := 0; i < mbody; i++ {
		ip := 6 * i
		ei := 0.5 * (math.Pow(yT[ip+4], 2) + math.Pow(yT[ip+5], 2) + math.Pow(yT[ip+6], 2))
		for j := i; j < mbody; j++ {
			jp := 6 * j
			dist := meps + math.Pow(yT[ip+1]-yT[jp+1], 2) + math.Pow(yT[ip+2]-yT[jp+2], 2) + math.Pow(yT[ip+3]-yT[jp+3], 2)
			ei = ei - m.mass[j]/math.Sqrt(dist)
		}
		en = en + m.mass[i]*ei
	}

	return en
}
