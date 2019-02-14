##  a bare bones Test-Suite for some critical functions
#using Base.Test
using Test, Random, LinearAlgebra, DelimitedFiles
include("unit_test_consts.jl")

@testset "All Tests" begin 

#Tests shrink/crossval other simpler helpers
@testset "Simple Helper Functions" begin
	@test true
end


@testset "POAP Tests" begin 
	#####Load up the POAP parameters
	include("../testHarness_Paper.jl")
	dat, header = readdlm("../Results/param_portExp_mtn2.csv", ',', header=true)
	Random.seed!(8675309)
	n = 2^10

	#Generate market
	thetas = vec(dat[1:n, 1])
	vs = vec(dat[1:n, 2])
	cs = vec(dat[1:n, 3])

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
	@test isapprox(thetaval, 0.15701062898375556) 	

	#FullInfo
	x_t[:] = x(cs, thetas)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.19235989290567318)

	#EB MLE
	tauMLE, x_t[:] = x_MLE(cs, muhat, vs)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.16258941230498414)
	@test isapprox(tauMLE,   0.3192500035088463)

	#EB MM
	tauMM, x_t[:] = x_MM(cs, muhat, vs)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.1632946586549207)
	@test isapprox(tauMM, 0.4235477938547216)

	#OR MSE
	x_t[:], tau_CV = x_OR_MSE(cs, muhat, thetas, vs)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.17011375672290682)
	@test isapprox(tau_CV, 1.1498678887082665)

	#SURE MSE
	x_t[:], tau_CV = x_sure_MSE(cs, muhat, vs)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.16884548740649977)
	@test isapprox(tau_CV, 0.8619343764734516)

	#Dirac Stein
	x_t[:], vals, objs = x_stein_exact(cs, muhat, vs, thetas)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.1693222241422801)
	for (ix, v) in enumerate(objs)
		@test isapprox(v, c_dirac_objs[ix], atol=1e-6)
	end

	#Box version of Stein
	h = n^-.16666
	x_t[:], vals, objs = x_stein_box(cs, muhat, vs, h, tau_step = .05)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.17053403283946592)
	for (ix, v) in enumerate(objs)
		@test isapprox(v, c_boxStein_objs[ix], atol=1e-6)
	end

	x_t[:], vals, objs = best_x_tau(cs, muhat, vs, thetas)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.17159684025607252)
	for (ix, v) in enumerate(objs)
		@test isapprox(v, c_oracleStein_objs[ix], atol=1e-6)
	end

	#Oracle Regs
	x_t[:], Gamma_grid, objs = KP.x_l2reg_CV(cs, muhat, vs, thetas, 
			Gamma_min=.1, Gamma_max=5, Gamma_step=.1)
	Gammahat = Gamma_grid[argmax(objs)]
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.17365203085789088)
	@test isapprox(Gammahat, 5.0)

	for (ix, v) in enumerate(objs)
		@test isapprox(v, c_oracleReg_objs[ix], atol=1e-6)
	end

	#Stein Reg
	x_t[:], Gamma_grid, objs = KP.x_stein_reg(cs, muhat, vs, 
											Gamma_min=.1, Gamma_max=.5, Gamma_step=.1)
	Gammahat = Gamma_grid[argmax(objs)]
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.16012121828686024)
	@test isapprox(Gammahat, .4)
	for (ix, v) in enumerate(objs)
		@test isapprox(v, c_SteinReg_objs[ix], atol=1e-6)
	end

	### Robust Methods
	thresh = sqrt(2*log(1/.1))
	x_t[:] = KP.x_robFW(cs, muhat, vs, thresh, TOL=1e-4)
	thetaval = dot(thetas, x_t)/n
	@test isapprox(thetaval, 0.15886998599283608)

	#LOO Validation
	x_t[:], Gamma_grid, objs = KP.x_LOO_reg(cs, muhat + noise, muhat - noise, vs, 
													Gamma_min=.1, Gamma_max=5, Gamma_step=.1)
	thetaval = dot(thetas, x_t)/n
	GammaLOO = Gamma_grid[argmax(objs)]

	@test isapprox(thetaval, 0.16415748804852356)
	@test isapprox(GammaLOO, 1.0)

	for (ix, v) in enumerate(objs)
		@test isapprox(v, c_LOOReg_objs[ix], atol=1e-6)
	end
end

#Simply tests the harness is runing without an errors
#Good to maintain signatures across functions
@testset "POAP Harness" begin
	include("../testHarness_Paper.jl")
    test_ReadData("../temp/temp_PortExp", 5, [100, 150], 8675309, "../Results/param_portExp_mtn2.csv")
    @test true  #to keep clean summary output
end

@testset "CLT Harness" begin
	include("../testCLTHarness.jl")
	test_POAPCLT("../Results/temp_POAPCLT", "../Results/param_portExp_mtn2.csv", 2, 100, [2 3], 8675309, "exponential", 1.0, 10., 1.)
	@test true #to keep clean summary output
end

end  #allTests stest set