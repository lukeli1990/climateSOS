function fLB(alpha)
    # alpha: 1 * (B*T);


    T=3

    B=14


#===========generate LB value =============#


mLB = Model(solver=GurobiSolver(OutputFlag=0))



@variable(mLB, 0<=aiwt[1:B,1:T]<=1)


@objective(mLB, Min, sum(alpha[i+B*(j-1)]*aiwt[i,j] for i=1:B,j=1:T))


status=solve(mLB)
#println("$status")
return getobjectivevalue(mLB)

end
