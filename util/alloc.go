package util

func MakeSquare(n uint) [][]float64 {
	return MakeRectangular(n, n)
}

func MakeRectangular(rows, cols uint) (rect [][]float64) {
	arr := make([]float64, rows*cols)
	rect = make([][]float64, rows)
	for i := range rect {
		rect[i] = arr[:cols]
		arr = arr[cols:]
	}
	return
}
