using JuMP, MathProgBase, Gurobi
m = Model(solver=GurobiSolver())
@variable(m, x[1:2,1:3])
@variable(m, y[1:2,1:3])
@constraint(m, x[1,2] + x[2,1] >= 1)
@constraint(m, x[1,2] + y[2,1] >= 1)
@constraint(m, x[1,3] + 2*y[1,3] <=3)
JuMP.build(m)
A=MathProgBase.getconstrmatrix(internalmodel(m))
