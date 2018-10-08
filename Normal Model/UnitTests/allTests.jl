##  a bare bones Test-Suite for some critical functions
using Base.Test
include("unit_test_consts.jl")

@testset "All Tests" begin 

#Tests shrink/crossval other simpler helpers
@testset "Simple Helper Functions" begin
	@test true
end;

@testset "POAP Tests" begin 
	#####Load up the POAP parameters
	include("../testHarness_Paper.jl")

	dat, header = readcsv("../Results/param_portExp_mtn1.csv", header=true)
	srand(8675309)
	const n = 2^10

	#Generate market of size n = 2^8
	thetas = dat[1:n, 1]
	vs = dat[1:n, 2]
	cs = dat[1:n, 3]
	cs /= quantile(cs, .2) #VG Remove this too on next fix.  Should be directly in the parameter.
	o = DefaultExp(cs, thetas, vs)
	muhat = zeros(Float64, n)
	noise = zeros(Float64, n)

	sim!(o, muhat)
	noise[:] = randn!(noise) ./ sqrt.(o.vs)

	###Now test each method.  Hope that testing objective is enough to catch most errors
	x_t = zeros(n)

	#SAA
	x_t[:] = x(cs, muhat)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.128596974016607) 	

	#FullInfo
	x_t[:] = x(cs, thetas)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.16429220417420648)

	#EB MLE
	tauMLE, x_t[:] = x_MLE(cs, muhat, vs)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.13522919274969386)
	@test isapprox(tauMLE,   0.3192500035088463)

	#EB MM
	tauMM, x_t[:] = x_MM(cs, muhat, vs)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.1358881373149644)
	@test isapprox(tauMM, 0.4235477938547216)

	#OR MSE
	x_t[:], tau_CV = x_OR_MSE(cs, muhat, thetas, vs)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.14222361529645114)
	@test isapprox(tau_CV, 1.1498678887082665)

	#SURE MSE
	x_t[:], tau_CV = x_sure_MSE(cs, muhat, vs)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.13999011139502265)
	@test isapprox(tau_CV, 0.8619343764734516)

	#Dirac Stein
	x_t[:], vals, objs = x_stein_exact(cs, muhat, vs, thetas)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.1419291050891903)
	for (ix, v) in enumerate(objs)
		@test isapprox(v, c_dirac_objs[ix], atol=1e-6)
	end

	#Box version of Stein
	const h = n^-.16666
	x_t[:], vals, objs = x_stein_box(cs, muhat, vs, h, tau_step = .05)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.1435738136345382)
	for (ix, v) in enumerate(objs)
		@test isapprox(v, c_boxStein_objs[ix], atol=1e-6)
	end

	x_t[:], vals, objs = best_x_tau(cs, muhat, vs, thetas)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.1438550512408222)
	for (ix, v) in enumerate(objs)
		@test isapprox(v, c_oracleStein_objs[ix], atol=1e-6)
	end

	#Oracle Regs
	x_t[:], Gamma_grid, objs = KP.x_l2reg_CV(cs, muhat, vs, thetas, 
			Gamma_min=.1, Gamma_max=5, Gamma_step=.1)
	Gammahat = Gamma_grid[indmax(objs)]
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.14620544329239127)
	@test isapprox(Gammahat, 5.0)

	for (ix, v) in enumerate(objs)
		@test isapprox(v, c_oracleReg_objs[ix], atol=1e-6)
	end

	#Stein Reg
	x_t[:], Gamma_grid, objs = KP.x_stein_reg(cs, muhat, vs, 
											Gamma_min=.1, Gamma_max=.5, Gamma_step=.1)
	Gammahat = Gamma_grid[indmax(objs)]
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.13190810832501446)
	@test isapprox(Gammahat, .4)

	for (ix, v) in enumerate(objs)
		@test isapprox(v, c_SteinReg_objs[ix], atol=1e-6)
	end

	### Robust Methods
	thresh = sqrt(2*log(1/.1))
	x_t[:] = KP.x_robFW(cs, muhat, vs, thresh, TOL=1e-4)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.1304026581578813)

	#LOO Validation
	x_t[:], Gamma_grid, objs = KP.x_LOO_reg(cs, muhat + noise, muhat - noise, vs, 
													Gamma_min=.1, Gamma_max=5, Gamma_step=.1)
	thetaval = dot(thetas, x_t)/n
	GammaLOO = Gamma_grid[indmax(objs)]
	@test isapprox(thetaval, 0.1374913242162149)
	@test isapprox(GammaLOO, 1.4)

	for (ix, v) in enumerate(objs)
		@test isapprox(v, c_LOOReg_objs[ix], atol=1e-6)
	end

end

end