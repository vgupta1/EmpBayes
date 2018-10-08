using Distributions, JuMP, Ipopt
include("KPNormal.jl")
using KP

function rob(cs, muhat, vs, r)
    const sqrt_vmin = sqrt(minimum(vs))
    function f(Gamma)
        xs, lam = KP.x_l2reg(cs, muhat, vs, Gamma)
#        println(Gamma, "\t", Gamma * sqrt_vmin * v_norm(xs, vs) - r)
        return Gamma * sqrt_vmin * v_norm(xs, vs) - r
    end
    println(f(.01), "\t", f(100))
    Gammastar = fzero(f, .01, 100)
    return KP.x_l2reg(cs, muhat, vs, Gammastar)
end

function rob_direct(cs, muhat, vs, r)
    m = Model(solver=IpoptSolver())
    @variable(m, xs[1:n] >= 0)
    @constraint(m, sum(cs[i]* xs[i] for i = 1:n) <= n)
    @NLobjective(m, Max, sum(muhat[i] * xs[i]/n for i = 1:n) - r/n * sqrt( sum(xs[i]^2 for i = 1:n)))
    solve(m)
end

v_norm(xs, vs) = sqrt(sum(xs.^2 ./ vs))

n =  1000
srand(8675309)
thetas = ones(n)
vs = ones(n)
cs = ones(n)

#From above
theta_l, v_l, c_l = 0.0, .01,  10.001
theta_m, v_m, c_m = 1.0, .75, 10.001
theta_h, v_h, c_h = 0.3,  2, 10.001

thetas[1:3:n] = theta_l
vs[1:3:n] = v_l
cs[1:3:n] = c_l

thetas[2:3:n] = theta_m
vs[2:3:n] = v_m
cs[2:3:n] = c_m

thetas[3:3:n] = theta_h
vs[3:3:n] = v_h
cs[3:3:n] = c_h

muhat = thetas + randn(n) ./ sqrt(vs)
sqrt_vmin = sqrt(minimum(vs))

tic()
xs, lam = rob(cs, muhat, vs, quantile(Normal(0, 1), .99))
toc()

tic()
xs2 = rob_direct(cs, muhat, vs, quantile(Normal(0, 1), .99)) 
toc()

