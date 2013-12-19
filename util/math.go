package util

import (
	"fmt"
	"math"
	"math/rand"
)

func Max(x, y int) int {
	if x > y {
		return x
	}
	return y
}

func Min(x, y int) int {
	if x < y {
		return x
	}
	return y
}

func ArrayEpsEquals(x, y []float64, eps float64) bool {
	if len(x) != len(y) {
		return false
	}
	for i := range x {
		if !EpsEqual(x[i], y[i], eps) {
			panic(fmt.Sprintln("Unequal entries at (%d): [%v, %v]", i, x[i], y[i]))
		}
	}
	return true
}

func EpsEqual(x, y, eps float64) bool {
	return math.Abs(x-y) < eps
}

func RandomInInterval(low, high float64) float64 {
	return low + (rand.Float64() * (high - low))
}
