# include("utils.jl")
using PyPlot
using Random
using LinearAlgebra
using Printf

# TODO: check GT solution given correct distribution
# using Distributions
# outliers = Laplace(mu_gt[1], 1)

rownorm(A) = mapslices(norm, A, dims=2)

p = rand(500,2) #< outliers

# --- create data
# Random.seed!(0)
p = vcat(p, .02*randn(200,2) .+ rand(1,2))
p = vcat(p, .02*randn(200,2) .+ rand(1,2))
p = vcat(p, .02*randn(200,2) .+ rand(1,2))
mu = rand(1,2)
mu_ls = sum(p, dims=1) ./ size(p,1)
r = ones(size(p,1),1)
pnorm = .2

mu_history = mu

for i=1:20
  global r, w, mu, mu_history
  
  # --- Compute residuals & weights
  r = rownorm(p.-mu)
  w = 1 ./ (r .+ 1e-3).^(2-pnorm)
  
  # --- Plotting weights
  clf()
  plot(mu[:,1], mu[:,2], "bx", markersize=7)
  plot(mu_history[:,1], mu_history[:,2], "b.-", markersize=4)
  plot(mu_ls[:,1], mu_ls[:,2], "ro", markersize=5)
  scatter(p[:,1], p[:,2], s=3, c=log.(w[:,1]), cmap=get_cmap("inferno"))
  gca().axis("square")
  gca().axis("off")
  xlim(0,1); ylim(0,1)
  display(gcf());

  filename = @sprintf "out/frame_%02.d" i
  savefig(filename, bbox_inches="tight")

  # --- Solve re-weighted LS problem to update the solution
  mu = sum(w .* p, dims=1) ./ sum(w, dims=1)
  mu_history = [mu_history; mu]
end

# --- Needs "brew install imagemagick"
run(`convert -delay 10 -loop 0 out/*.png animation.gif`)

# --- To merge them all into one (at end)
# run(`convert g*.gif out.gif`)