/**
 * Created with IntelliJ IDEA.
 * User: Georg
 * Date: 21.10.13
 * Time: 13:31
 * To change this template use File | Settings | File Templates.
 */
package problems

type Problem interface {
	Description() string
	Initialize() []float64
	Fcn(t float64, yT []float64, dy_out []float64)
}

type TiledProblem interface {
	Problem
	FcnBlock(startIdx, blockSize uint, t float64, yT []float64, dy_out []float64)
}
