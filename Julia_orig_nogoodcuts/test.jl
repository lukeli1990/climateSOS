using JuMP, MathProgBase, Gurobi
m = Model(solver=GurobiSolver())
@variable(m, x[1:2])
@constraint(m, x[1] + x[2] >= 1)
@constraint(m, x[1] + 2*x[2] <=3)
JuMP.build(m)
MathProgBase.getconstrmatrix(internalmodel(m))

#include("first_stage.jl")
#x0=1
#costval=[1 -1 3]
#x=f1st(x0,costval)
