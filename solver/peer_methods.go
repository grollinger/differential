package solver

import (
	"errors"
	"math"
)

type PeerMethod uint

const (
	EPP2 = PeerMethod(iota)
	EPP4
	EPP4y2
	EPP4y3
	EPP4_06809
	EPP6p1
	EPP6j1
	EPP8_d
	EPP8sp8
	EPP_x1
	EPP_x2
	NumberOfPeerMethods = uint(iota)
)

func NewPeer(m PeerMethod) (i Integrator, err error) {
	var p peer
	p.method = m
	err = p.setCoeffs()
	if err != nil {
		p = peer{}
	}
	i = &p

	return
}

//-- solves Vandermonde systems PM_new*V=PM
func vanderMonde(pc []float64, pm [][]float64) {
	var i, j, k int
	n := len(pc)

	for k = 0; k < n-1; k++ {
		for j = n - 1; j >= k+1; j-- {
			for i = 0; i < n; i++ {
				pm[i][j] = pm[i][j] - pm[i][j-1]*pc[k]
			}
		}
	}

	for k = n - 2; k >= 0; k-- {
		for j = k + 1; j < n; j++ {
			for i = 0; i < n; i++ {
				pm[i][j] = pm[i][j] / (pc[j] - pc[j-k-1])
			}
		}

		for j = k; j < n-1; j++ {
			for i = 0; i < n; i++ {
				pm[i][j] = pm[i][j] - pm[i][j+1]
			}
		}
	}
}

func (p *peer) setCoeffs() (err error) {
	switch p.method {
	case EPP2:
		p.setEPP2Coeffs()
	case EPP4:
		p.setEPP4Coeffs()
	case EPP4y2:
		p.order, p.stages, p.sigx = 4, 4, 1.6
		p.name = "EPP4y2"
		p.allocateCoeffs()
		p.setEPP4y2Coeffs()
	case EPP4y3:
		p.order, p.stages, p.sigx = 4, 4, 1.6
		p.ema = 0.3125
		p.name = "EPPy3"
		p.allocateCoeffs()
		p.setEPP4y3Coeffs()
	case EPP4_06809:
		p.order, p.stages, p.sigx = 4, 4, 1.6
		p.name = "EPP4_06809"
		p.allocateCoeffs()
		p.setEPP4_06809Coeffs()
	case EPP6p1:
		p.order, p.stages, p.sigx = 6, 6, 1.5
		p.ema = 0.125
		p.name = "EPP6p1"
		p.allocateCoeffs()
		p.setEPP6p1Coeffs()
	case EPP6j1:
		p.order, p.stages, p.sigx = 6, 6, 1.5
		p.ema = 0.125
		p.name = "EPP6j1"
		p.allocateCoeffs()
		p.setEPP6j1Coeffs()
	case EPP8_d:
		p.order, p.stages, p.sigx = 8, 8, 1.4
		p.name = "EPP8_d"
		p.allocateCoeffs()
		p.setEPP8_dCoeffs()
	case EPP8sp8:
		p.order, p.stages, p.sigx = 8, 8, 1.4
		p.name = "EPP8sp8"
		p.allocateCoeffs()
		p.setEPP8sp8Coeffs()
	case EPP_x1: // absc=-0.524, B-norm=0.35*8
		p.order, p.stages, p.sigx = 8, 8, 1.4
		p.name = "EPP_x1"
		p.allocateCoeffs()
		p.setEPP_x1Coeffs()
	case EPP_x2: // absc=-0.466, efc=0, B-norm=0.35*8 /20061101
		p.order, p.stages, p.sigx = 8, 8, 1.4
		p.ema = 1.0 // interval [0.2,2]
		p.name = "EPP_x2"
		p.allocateCoeffs()
		p.setEPP_x2Coeffs()
	default:
		err = errors.New("unknown peer method")
		return
	}

	p.findMinMaxNodes()

	p.ensureOneRowSums()

	// auxiliary parameters for local error model
	p.ema = 0.0 //ausschalten
	p.emaoh = math.Pow(p.ema, float64(p.order)/2.0)

	var i, j, k uint

	//compute matrix PCV = diag(C)*VdM
	for i = 0; i < p.stages; i++ {
		p.cv[i][0] = p.c[i]
		for j = 1; j < p.stages; j++ {
			p.cv[i][j] = p.c[i] * p.cv[i][j-1]
		}
	}

	// compute matrix PA0, ppv used as temporary memory
	for i = 0; i < p.stages; i++ {
		for j = 0; j < p.stages; j++ {
			p.pv[i][j] = -p.b[i][j]
		}
	}

	// now \ones*e_s^T-B
	for i = 0; i < p.stages; i++ {
		p.pv[i][p.stages-1] = p.pv[i][p.stages-1] + 1.0
	}

	//pa0 = matmul(ppv,pcv)
	for i = 0; i < p.stages; i++ {
		for j = 0; j < p.stages; j++ {
			// unnecessary because of go memory initializer
			// p.a0[i][j] = 0.0
			for k = 0; k < p.stages; k++ {
				p.a0[i][j] = p.a0[i][j] + p.pv[i][k]*p.cv[k][j]
			}
		}
	}

	// scale columns
	for j = 1; j < p.stages; j++ {
		for i = 0; i < p.stages; i++ {
			p.a0[i][j] = p.a0[i][j] / (float64(j) + 1.0)
		}
	}

	vanderMonde(p.c, p.a0) // now PA0=(\ones*e_s^T-B)*C*V*(V*D)^(-1)

	// compute matrix PPV:
	for j = 0; j < p.stages; j++ {
		p.pv[0][j] = 1.0
	}

	for i = 1; i < p.stages; i++ {
		p.pv[i][0] = 0.0
		for j = 0; j < p.stages-1; j++ {
			p.pv[i][j+1] = p.pv[i][j] + p.pv[i-1][j]
		}
	}

	// scale rows
	for i = 0; i < p.stages; i++ {
		for j = 0; j < p.stages; j++ {
			p.pv[i][j] = p.pv[i][j] / (float64(i) + 1.0)
		}
	}

	vanderMonde(p.c, p.pv) // now PPV=D^(-1)*P*V^(-1)

	// error estimate with last row of PPV
	copy(p.e, p.pv[p.stages-1])

	return
}

func stagesOf(m PeerMethod) uint {
	return uint(m) % 10
}

func (p *peer) ensureOneRowSums() {
	// row sums of B must be 1 exactly
	var i, j uint
	for i = 0; i < p.stages; i++ {
		s := 1.0
		for j = 0; j < p.stages; j++ {
			s = s - p.b[i][j]
		}
		p.b[i][p.stages-1] = p.b[i][p.stages-1] + s
	}
}

func (p *peer) findMinMaxNodes() {
	// ! minimal and maximal nodes:
	p.icmin = 0
	p.icmax = uint(p.stages) - 1

	var i uint
	for i = 0; i < p.stages; i++ {
		if p.c[i] < p.c[p.icmin] {
			p.icmin = i
		}
		if p.c[i] > p.c[p.icmax] {
			p.icmax = i
		}
	}
}

func (p *peer) allocateCoeffs() {
	p.c, p.e = make([]float64, p.stages), make([]float64, p.stages)
	p.b = makeSquare(p.stages)
	p.a0 = makeSquare(p.stages)
	p.cv = makeSquare(p.stages)
	p.pv = makeSquare(p.stages)
}

func (p *peer) setEPP2Coeffs() {
	// Fortran Code says order = 4 ... really?
	p.order, p.stages, p.sigx = 2, 2, 1.5
	p.name = "EPP2"
	p.allocateCoeffs()
	// p.ema = 0.0

	p.c[0] = -1.0
	p.c[1] = 1.0

	p.b[0][0], p.b[0][1] = 0.5, 0.5
	p.b[1][0], p.b[1][1] = 0.5, 0.5
}

func (p *peer) setEPP4Coeffs() {
	p.order, p.stages, p.sigx = 4, 4, 1.4
	p.name = "EPP4"
	p.allocateCoeffs()
	// p.ema = 0.0

	p.c[0] = -1.0
	p.c[1] = -2.0 / 5.0
	p.c[2] = 11.0 / 20.0
	p.c[3] = 1.0

	p.b[0][0] = -0.810068472
	p.b[0][1] = -0.437909571
	p.b[0][2] = 5.618482265
	p.b[0][3] = -3.370504222
	p.b[1][0] = 0.008525097
	p.b[1][1] = -1.665642634
	p.b[1][2] = 5.799902598
	p.b[1][3] = -3.142785061
	p.b[2][0] = -0.476618420
	p.b[2][1] = 0.5422323993
	p.b[2][2] = 1.087151458
	p.b[2][3] = -0.1527654373
	p.b[3][0] = -1.468520018
	p.b[3][1] = 3.36438373
	p.b[3][2] = -3.284423359
	p.b[3][3] = 2.388559647
}

func (p *peer) setEPP4y2Coeffs() {
	// p.ema = 0.0

	p.c[0] = 0.44856672599000208
	p.c[1] = 1.39573694851711427
	p.c[2] = 1.862002092540212
	p.c[3] = 1.0

	p.b[0][0] = -0.00085592109218945
	p.b[0][1] = -0.00002375632236906
	p.b[0][2] = 0.00011240281897356
	p.b[0][3] = 1.00076727459558495
	p.b[1][0] = -0.18463151693356136
	p.b[1][1] = 0.0166248564168753
	p.b[1][2] = -0.08694965174179996
	p.b[1][3] = 1.25495631225848602
	p.b[2][0] = -0.01340990521951441
	p.b[2][1] = 0.00305657050832219
	p.b[2][2] = -0.01576893532468585
	p.b[2][3] = 1.02612227003587808
	p.b[3][0] = 0.0
	p.b[3][1] = 0.0
	p.b[3][2] = 0.0
	p.b[3][3] = 1.0
}

func (p *peer) setEPP4y3Coeffs() {

	p.c[0] = 1.33880820864483004
	p.c[1] = 1.70380840062134099
	p.c[2] = 1.86823097147835395
	p.c[3] = 1.0

	p.b[0][0] = 0.00534248151938605
	p.b[0][1] = -0.12787556805087963
	p.b[0][2] = -0.01346305501819259
	p.b[0][3] = 1.13599614154968617
	p.b[1][0] = 0.00291846709251725
	p.b[1][1] = -0.00273308699096256
	p.b[1][2] = -0.07182299833148626
	p.b[1][3] = 1.07163761822993158
	p.b[2][0] = 0.00028874763346411
	p.b[2][1] = -0.00495214157424223
	p.b[2][2] = -0.00260939452842348
	p.b[2][3] = 1.00727278846920161
	p.b[3][0] = 0.0
	p.b[3][1] = 0.0
	p.b[3][2] = 0.0
	p.b[3][3] = 1.0
}
func (p *peer) setEPP4_06809Coeffs() {

	p.c[0] = -1.067193866512852
	p.c[1] = -2.756684444690223e-1
	p.c[2] = 1.946690102724974
	p.c[3] = 1.0

	p.b[0][0] = -2.716560064574534e-1
	p.b[0][1] = 8.891878323750345e-1
	p.b[0][2] = -2.925373526104706e-1
	p.b[0][3] = 6.750055266928895e-1
	p.b[1][0] = 1.349945068349794e-1
	p.b[1][1] = 2.852936122417532e-1
	p.b[1][2] = -2.183244224563100e-1
	p.b[1][3] = 7.980363033795774e-1
	p.b[2][0] = -1.316829664799590e-1
	p.b[2][1] = 1.684954092016832
	p.b[2][2] = 1.085405981287241
	p.b[2][3] = -1.638677106824114
	p.b[3][0] = 2.358621659008970e-1
	p.b[3][1] = 5.306425158463548e-1
	p.b[3][2] = 3.325389053242893e-1
	p.b[3][3] = -9.904358707154119e-2
}

func (p *peer) setEPP6p1Coeffs() {

	p.c[0] = -1.31059599683912621
	p.c[1] = 1.97665537290660046
	p.c[2] = 1.6649784136534037
	p.c[3] = 1.1586933567385133
	p.c[4] = 0.5741070138915608
	p.c[5] = 1.0

	p.b[0][0] = 0.13895408777711341
	p.b[0][1] = -0.08837957588373672
	p.b[0][2] = 0.32222771745172996
	p.b[0][3] = -0.33501322375446315
	p.b[0][4] = 0.10714669158143159
	p.b[0][5] = 0.8550643028279249
	p.b[1][0] = 0.76848393019728558
	p.b[1][1] = -0.39601237423788838
	p.b[1][2] = 1.71331670705600939
	p.b[1][3] = -0.16709652789611371
	p.b[1][4] = 0.54842536517292088
	p.b[1][5] = -1.46711710029221376
	p.b[2][0] = 0.17348417746015228
	p.b[2][1] = -0.07154407321571946
	p.b[2][2] = 0.29463725844330063
	p.b[2][3] = -0.08003079966892302
	p.b[2][4] = -0.21195447314566
	p.b[2][5] = 0.89540791012684957
	p.b[3][0] = 0.00250303712641231
	p.b[3][1] = 0.00275619078550659
	p.b[3][2] = -0.02788674768628258
	p.b[3][3] = -0.06614667221305652
	p.b[3][4] = -0.12645307545127383
	p.b[3][5] = 1.21522726743869403
	p.b[4][0] = 0.00063646601897938
	p.b[4][1] = -0.00344178637599682
	p.b[4][2] = 0.01015487188250037
	p.b[4][3] = -0.02819421169069512
	p.b[4][4] = 0.02856770023053085
	p.b[4][5] = 0.99227695993468134
	p.b[5][0] = 0.0
	p.b[5][1] = 0.0
	p.b[5][2] = 0.0
	p.b[5][3] = 0.0
	p.b[5][4] = 0.0
	p.b[5][5] = 1.0
}
func (p *peer) setEPP6j1Coeffs() {

	p.c[0] = 6.1182488158460324e-1
	p.c[1] = 1.0734784354567433
	p.c[2] = 1.7733348046756701
	p.c[3] = 1.9723174701317718
	p.c[4] = 1.4155260278449762
	p.c[5] = 1.0

	p.b[0][0] = -2.0180146181687607e-4
	p.b[0][1] = 1.6304614802106130e-2
	p.b[0][2] = -1.2851544816318155e-2
	p.b[0][3] = 2.6210658256845634e-3
	p.b[0][4] = 3.7912816867901877e-3
	p.b[0][5] = 9.9033638396355417e-1
	p.b[1][0] = 5.4402447198783130e-5
	p.b[1][1] = 2.4618095745274586e-4
	p.b[1][2] = 4.1283950478615055e-3
	p.b[1][3] = 1.0819600004066808e-3
	p.b[1][4] = -6.4960842108610951e-3
	p.b[1][5] = 1.0009851457579413
	p.b[2][0] = -1.3283089416477777e-4
	p.b[2][1] = 2.6585250025434079e-4
	p.b[2][2] = -3.9985754350926959e-4
	p.b[2][3] = -1.4512932179574525e-2
	p.b[2][4] = 1.0221144035648493e-2
	p.b[2][5] = 1.0045586240813458
	p.b[3][0] = 1.9579814800354238e-4
	p.b[3][1] = -1.4979121212198304e-4
	p.b[3][2] = 1.4147303648950299e-4
	p.b[3][3] = -1.8142965943514442e-4
	p.b[3][4] = -1.7806845742653309e-2
	p.b[3][5] = 1.0178007954297175
	p.b[4][0] = -7.6319822223659301e-6
	p.b[4][1] = 1.8177963233112666e-4
	p.b[4][2] = -1.7553812394822423e-4
	p.b[4][3] = -4.0614157234655195e-5
	p.b[4][4] = 5.3690770730854430e-4
	p.b[4][5] = 9.9950509692376555e-1
	p.b[5][0] = 0.0
	p.b[5][1] = 0.0
	p.b[5][2] = 0.0
	p.b[5][3] = 0.0
	p.b[5][4] = 0.0
	p.b[5][5] = 1.0
}
func (p *peer) setEPP8_dCoeffs() {
	p.c[0] = 0.26041740957753135
	p.c[1] = 0.52923626823623069
	p.c[2] = 1.54653689839871537
	p.c[3] = 1.77514379674033294
	p.c[4] = 0.090107569010498
	p.c[5] = 0.19995116126552217
	p.c[6] = 1.32615523512850911
	p.c[7] = 1.0

	p.b[0][0] = 0.00218540643609441
	p.b[0][1] = 0.00964151773849017
	p.b[0][2] = 0.02740905027639759
	p.b[0][3] = -0.00757481375423608
	p.b[0][4] = -0.00627931764830473
	p.b[0][5] = -0.00651613228833201
	p.b[0][6] = -0.00103882981570947
	p.b[0][7] = 0.98217311905560012
	p.b[1][0] = 0.0026911855247508
	p.b[1][1] = -0.00999004293842317
	p.b[1][2] = 0.01299877362586784
	p.b[1][3] = -0.00362117297559769
	p.b[1][4] = -0.01393348294082885
	p.b[1][5] = 0.00206411701387576
	p.b[1][6] = -0.00794032918151776
	p.b[1][7] = 1.01773095187187307
	p.b[2][0] = 0.00199608410884834
	p.b[2][1] = -0.01002762418952864
	p.b[2][2] = 0.02034106943279384
	p.b[2][3] = -0.00932460631409766
	p.b[2][4] = -0.00937443352246159
	p.b[2][5] = 0.00649579387949097
	p.b[2][6] = -0.02592014846435581
	p.b[2][7] = 1.02581386506931055
	p.b[3][0] = 0.00103036221311367
	p.b[3][1] = -0.01190896158535306
	p.b[3][2] = 0.02063699940945861
	p.b[3][3] = -0.00672959115904196
	p.b[3][4] = -0.00583850891833535
	p.b[3][5] = 0.03351758481419382
	p.b[3][6] = -0.03846471330760932
	p.b[3][7] = 1.0077568285335736
	p.b[4][0] = 0.00148756342489522
	p.b[4][1] = -0.01047547599327585
	p.b[4][2] = 0.0208838878315454
	p.b[4][3] = -0.00663219960933854
	p.b[4][4] = -0.00563964345840616
	p.b[4][5] = 0.02716167598460901
	p.b[4][6] = -0.00664272798306554
	p.b[4][7] = 0.97985691980303647
	p.b[5][0] = 0.00205269876381091
	p.b[5][1] = -0.00960860565429818
	p.b[5][2] = 0.0205427263285296
	p.b[5][3] = -0.00625381166734364
	p.b[5][4] = -0.00577560681234781
	p.b[5][5] = 0.01293837851580094
	p.b[5][6] = 0.00660793064341287
	p.b[5][7] = 0.97949628988243531
	p.b[6][0] = 0.00184988219096286
	p.b[6][1] = -0.01079724930569055
	p.b[6][2] = 0.02051636130170474
	p.b[6][3] = -0.00632847064898024
	p.b[6][4] = -0.00550013664675566
	p.b[6][5] = 0.01390490619566548
	p.b[6][6] = -0.01323688831455224
	p.b[6][7] = 0.99959159522764561
	p.b[7][0] = 0.00186731024859287
	p.b[7][1] = -0.01061165690522434
	p.b[7][2] = 0.02042910717986102
	p.b[7][3] = -0.00630049368165016
	p.b[7][4] = -0.00568571957918147
	p.b[7][5] = 0.01348102696866063
	p.b[7][6] = -0.01331088571679290
	p.b[7][7] = 1.00013131148573434
}
func (p *peer) setEPP8sp8Coeffs() {
	p.c[0] = 0.70541387781778147
	p.c[1] = 1.30641486071640562
	p.c[2] = 0.31760983680370311
	p.c[3] = 1.34289163028271588
	p.c[4] = 0.26801849628336664
	p.c[5] = 1.06975280780012286
	p.c[6] = 1.20363542351617524
	p.c[7] = 1.0

	p.b[0][0] = 0.27218496097888581
	p.b[0][1] = 0.24481602761226617
	p.b[0][2] = 0.14324265273812842
	p.b[0][3] = 0.07518091759856021
	p.b[0][4] = 0.10880381601605394
	p.b[0][5] = -0.08559047504042014
	p.b[0][6] = -0.89848704629310313
	p.b[0][7] = 1.1398491463896287
	p.b[1][0] = 0.27218494251985784
	p.b[1][1] = 0.24475826229310614
	p.b[1][2] = 0.14312713840849838
	p.b[1][3] = 0.07511662694065275
	p.b[1][4] = 0.10880496489216468
	p.b[1][5] = -0.08573012564915534
	p.b[1][6] = -0.89863712500172939
	p.b[1][7] = 1.14037531559660494
	p.b[2][0] = 0.27218494691659497
	p.b[2][1] = 0.24475827494648921
	p.b[2][2] = 0.14324960884617646
	p.b[2][3] = 0.07514374234973973
	p.b[2][4] = 0.10870161200141529
	p.b[2][5] = -0.08562615141184479
	p.b[2][6] = -0.89856217679317559
	p.b[2][7] = 1.14015014314460471
	p.b[3][0] = 0.2721849293851366
	p.b[3][1] = 0.24475827854943941
	p.b[3][2] = 0.1432496044378971
	p.b[3][3] = 0.07516501652473434
	p.b[3][4] = 0.1086230286174295
	p.b[3][5] = -0.08576394151117088
	p.b[3][6] = -0.89869867862268446
	p.b[3][7] = 1.1404817626192184
	p.b[4][0] = 0.27218493682276326
	p.b[4][1] = 0.2447582729025796
	p.b[4][2] = 0.14324957805862669
	p.b[4][3] = 0.07516500809749052
	p.b[4][4] = 0.10873535122317697
	p.b[4][5] = -0.08577918038240442
	p.b[4][6] = -0.89851969781206312
	p.b[4][7] = 1.1402057310898305
	p.b[5][0] = 0.27218494476634404
	p.b[5][1] = 0.24475826960211763
	p.b[5][2] = 0.14324957747511375
	p.b[5][3] = 0.07516500706866706
	p.b[5][4] = 0.10873538294098683
	p.b[5][5] = -0.08567964439138113
	p.b[5][6] = -0.89843186740710734
	p.b[5][7] = 1.14001832994525916
	p.b[6][0] = 0.27218494660576334
	p.b[6][1] = 0.2447582805123369
	p.b[6][2] = 0.14324957501518646
	p.b[6][3] = 0.0751650064576131
	p.b[6][4] = 0.10873538290916957
	p.b[6][5] = -0.08567965711914102
	p.b[6][6] = -0.89853485595644752
	p.b[6][7] = 1.14012132157551916
	p.b[7][0] = 0.27218494660493086
	p.b[7][1] = 0.24475827948013657
	p.b[7][2] = 0.14324958032207845
	p.b[7][3] = 0.07516500887490214
	p.b[7][4] = 0.10873538145127106
	p.b[7][5] = -0.08567966433201819
	p.b[7][6] = -0.89853483288304983
	p.b[7][7] = 1.14012130048174894
}
func (p *peer) setEPP_x1Coeffs() {
	p.c[0] = -1.020253410235809
	p.c[1] = -7.973369854084624e-1
	p.c[2] = -5.523527869930042e-1
	p.c[3] = -1.601795297055826e-1
	p.c[4] = 2.567353840060831e-1
	p.c[5] = 5.758479688560444e-1
	p.c[6] = 8.048768575883170e-1
	p.c[7] = 1.0

	p.b[0][0] = 2.625189510221092e-5
	p.b[0][1] = -1.898215388937858e-5
	p.b[0][2] = 1.817009753212212e-3
	p.b[0][3] = 9.967928831100421e-1
	p.b[0][4] = 1.741824921921685e-4
	p.b[0][5] = 2.174441508324741e-3
	p.b[0][6] = -1.677221467830472e-4
	p.b[0][7] = -7.980644582009861e-4
	p.b[1][0] = 8.100101125281656e-4
	p.b[1][1] = 1.002801199171737e-3
	p.b[1][2] = -1.118614707367094e-4
	p.b[1][3] = 9.970279280278038e-1
	p.b[1][4] = -1.983509579940047e-3
	p.b[1][5] = 7.804032379917670e-4
	p.b[1][6] = 7.881153868976325e-4
	p.b[1][7] = 1.686113086283548e-3

	p.b[2][0] = -9.879713512756130e-4
	p.b[2][1] = 8.020695098409326e-5
	p.b[2][2] = 1.666760283706175e-3
	p.b[2][3] = 9.953316156454997e-1
	p.b[2][4] = 2.508073453472758e-4
	p.b[2][5] = -3.654635962011422e-4
	p.b[2][6] = 1.959602069475001e-3
	p.b[2][7] = 2.064442652464573e-3
	p.b[3][0] = -8.835912374043340e-4
	p.b[3][1] = -1.760815796251268e-3
	p.b[3][2] = 2.331454645143357e-3
	p.b[3][3] = 9.981121059461671e-1
	p.b[3][4] = 2.273217656663430e-3
	p.b[3][5] = -6.558878934506733e-5
	p.b[3][6] = 1.101322224606321e-4
	p.b[3][7] = -1.169146474339194e-4
	p.b[4][0] = 2.214901854024395e-3
	p.b[4][1] = -4.929384902160960e-4
	p.b[4][2] = 1.412225706887807e-3
	p.b[4][3] = 9.961893228068536e-1
	p.b[4][4] = -1.921087072457036e-4
	p.b[4][5] = -2.164510616598476e-3
	p.b[4][6] = 1.669222622278326e-3
	p.b[4][7] = 1.363884824016043e-3
	p.b[5][0] = 1.576721961555554e-3
	p.b[5][1] = 1.317624066616264e-4
	p.b[5][2] = 1.386449758969554e-3
	p.b[5][3] = 9.934659543403941e-1
	p.b[5][4] = 1.424864416795372e-3
	p.b[5][5] = 9.666678111011650e-4
	p.b[5][6] = -1.288319653300130e-3
	p.b[5][7] = 2.335898957822821e-3
	p.b[6][0] = -1.429351554370223e-3
	p.b[6][1] = -2.054822802735597e-4
	p.b[6][2] = 5.036279143570400e-4
	p.b[6][3] = 1.001130852115276
	p.b[6][4] = -7.867597108633510e-4
	p.b[6][5] = 1.592378375989060e-3
	p.b[6][6] = 2.548549580728228e-4
	p.b[6][7] = -1.060119818187650e-3
	p.b[7][0] = 6.277247319384195e-4
	p.b[7][1] = 4.892457115011240e-4
	p.b[7][2] = -1.541906260922679e-3
	p.b[7][3] = 1.001178331446194
	p.b[7][4] = 1.541403411299333e-3
	p.b[7][5] = -1.715851337805785e-3
	p.b[7][6] = -1.947535937192549e-3
	p.b[7][7] = 1.368588234988504e-3
}
func (p *peer) setEPP_x2Coeffs() {

	p.c[0] = -1.514542417302030
	p.c[1] = -1.003995798476134
	p.c[2] = -5.372648489667148e-1
	p.c[3] = -6.969892603966277e-2
	p.c[4] = 2.867223168569619e-1
	p.c[5] = 7.943872874460213e-1
	p.c[6] = 1.516388473646935
	p.c[7] = 1.0

	p.b[0][0] = -2.806877716477130e-2
	p.b[0][1] = 1.135981190423357e-2
	p.b[0][2] = -2.110383374673507e-2
	p.b[0][3] = 1.037885271481855
	p.b[0][4] = -2.490268975167300e-2
	p.b[0][5] = 7.668414990349029e-3
	p.b[0][6] = 1.814611059875396e-2
	p.b[0][7] = -9.843083120119264e-4
	p.b[1][0] = 3.833268630087496e-3
	p.b[1][1] = -2.144635474502071e-2
	p.b[1][2] = -3.268424349441350e-2
	p.b[1][3] = 1.061028463256470
	p.b[1][4] = -1.543617522327160e-2
	p.b[1][5] = 2.390602005607571e-2
	p.b[1][6] = -4.431530991242315e-3
	p.b[1][7] = -1.476944748868486e-2
	p.b[2][0] = 1.517768428061389e-2
	p.b[2][1] = -4.118587370592280e-3
	p.b[2][2] = -3.428305617434239e-3
	p.b[2][3] = 9.841348202355257e-1
	p.b[2][4] = -2.899959231656522e-2
	p.b[2][5] = 3.010638273007045e-2
	p.b[2][6] = 2.187082102577848e-2
	p.b[2][7] = -1.474322296739694e-2
	p.b[3][0] = 2.774930716975024e-2
	p.b[3][1] = 3.233760768362689e-2
	p.b[3][2] = 1.578243438357346e-2
	p.b[3][3] = 8.980551794364540e-1
	p.b[3][4] = 1.081300247363562e-3
	p.b[3][5] = 1.959734287518365e-2
	p.b[3][6] = -1.255493065641702e-3
	p.b[3][7] = 6.652321269689837e-3
	p.b[4][0] = 6.882004136347424e-3
	p.b[4][1] = -2.207359291833351e-2
	p.b[4][2] = 1.806414842816704e-2
	p.b[4][3] = 9.890131245336967e-1
	p.b[4][4] = 5.357020858953068e-3
	p.b[4][5] = -9.626265478966035e-3
	p.b[4][6] = -1.809330236593781e-2
	p.b[4][7] = 3.047686280607328e-2
	p.b[5][0] = -1.803090133918039e-2
	p.b[5][1] = -4.588805985737581e-2
	p.b[5][2] = -2.519561283870134e-2
	p.b[5][3] = 1.075261749109603
	p.b[5][4] = 2.031592200888563e-2
	p.b[5][5] = -3.455936510468920e-2
	p.b[5][6] = 6.132333433229704e-3
	p.b[5][7] = 2.196393458822805e-2
	p.b[6][0] = 4.751715502768429e-3
	p.b[6][1] = -1.080836510327213e-2
	p.b[6][2] = 4.242786141675581e-3
	p.b[6][3] = 9.614609736912729e-1
	p.b[6][4] = 2.758883811192216e-2
	p.b[6][5] = 7.255102679141721e-3
	p.b[6][6] = -9.507573358168125e-3
	p.b[6][7] = 1.501652233465931e-2
	p.b[7][0] = 1.700680076802063e-2
	p.b[7][1] = 8.759888270257397e-3
	p.b[7][2] = -1.497101816869549e-2
	p.b[7][3] = 9.355634245613111e-1
	p.b[7][4] = 1.410966089361167e-2
	p.b[7][5] = 2.494725465581370e-2
	p.b[7][6] = 1.866148320743530e-2
	p.b[7][7] = -4.077494187754325e-3
}
