function f1st(x0,costval,UB,scenariow,loadbus,loadnum,ei,di,iter,oldaiwt,orderscenw)
    # x0 UB iter in B*1; costval in B*T;
    # CANT DECOMPOSE and Need cut in each iteration
    # input loadbus will determine a relationship between fw and aiwt

    #system parameter
    T=3 #time interval
    B=14 #bus number
    sizew=20 #size of scenario set

m1st = Model(solver=GurobiSolver(OutputFlag=0))


@variable(m1st, x0[i]<=x[i=1:B,1:T]<= UB[i],Int)
@variable(m1st, fw[1:sizew],Bin)
@variable(m1st, aiwt[1:B,1:T,1:sizew],Bin)
@variable(m1st, 0<=aiwtxit[i=1:B,1:T,1:sizew]<=UB[i]) #slack variable for product of awit and xit

@objective(m1st, Min, sum(costval[i,j]*x[i,j] for i=1:B,j=1:T))

#initial condition = x0
@constraint(m1st, [i=1:B], x[i,1]==x0[i])

#increasing adaptation at each time step
@constraint(m1st, [i=1:B,t=2:T], x[i,t]>=x[i,t-1])

#sum of total fw
@constraint(m1st, sum(fw[i] for i=1:sizew)>=0.8*sizew) #90% 60% 20% (difficult to solve)

#link between fw and a
for i=1:loadnum

    @constraint(m1st, [j=1:T, w=1:sizew], aiwt[loadbus[i],j,w]>=fw[w])

end

#constraint 13a with slack
@constraint(m1st, [i=1:B,j=1:T,w=1:sizew], (2*aiwt[i,j,w]-1)*scenariow[i,j,w] + ei[i]+di[i]*x[i,j] - 2*di[i]*aiwtxit[i,j,w] <=0)


#mccormick for aiwtxit
for i=1:B
    for j=1:T
        for w=1:sizew
            mccormick(m1st,aiwtxit[i,j,w],x[i,j],aiwt[i,j,w],x0[i],UB[i],0,1)
        end
    end
end

#no good cut
if iter>1
    for indcut=1:iter-1
        @constraint(m1st,sum((1-oldaiwt[i,j,w,indcut])*aiwt[i,j,w] + (1-aiwt[i,j,w])*oldaiwt[i,j,w,indcut] for i=1:B,j=1:T,w=1:sizew)>=1 )
    end
end

#order on aiwt
@constraint(m1st,[i=1:B,t=1:T],aiwt[i,t,orderscenw[i,t,sizew]]-aiwt[i,t,orderscenw[i,t,1]]<=1)

@constraint(m1st,[i=1:B,t=1:T,w=2:sizew],aiwt[i,t,orderscenw[i,t,w]]>= aiwt[i,t,orderscenw[i,t,w-1]])



# Remark: 14 (g) is not needed anymore
status = solve(m1st)
println("$status")
xsol=getvalue(x) # column T*1
fwsol=getvalue(fw)
aiwtsol=getvalue(aiwt)
return xsol,fwsol,aiwtsol
end
