function f1st(x0,costval,UB,scenariow,loadbus,loadnum,ei,di,iter,count_fw_vio,alpha_col,fw_vio_ind,hvalue,LBvalue,conflvl)
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
@constraint(m1st, sum(fw[i] for i=1:sizew)>=sizew*conflvl) #90 percent requirement 60 20 40
# @constraint(m1st,fw[19]==1)
# @constraint(m1st,fw[14]==1)
# @constraint(m1st,fw[13]==1)
# @constraint(m1st,fw[9]==1)

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

if iter>1
    @constraint(m1st,[j=1:(iter-1),i=1:count_fw_vio[j]], alpha_col[:,i,j]'*[aiwt[:,1,fw_vio_ind[i,j]];aiwt[:,2,fw_vio_ind[i,j]];aiwt[:,3,fw_vio_ind[i,j]]] + (hvalue[i,j]-LBvalue[i,j])*(1-fw[fw_vio_ind[i,j]])>=hvalue[i,j])
    #@constraint(m1st,[i=1:count_fw_vio[iter-1],j=1:(iter-1)], alpha_col[:,i,j]'*[x[:,1];x[:,2];x[:,3]]>=cutcons[i,j]-1e-3 )
end




# Remark: 14 (g) is not needed anymore
status = solve(m1st)
println("$status")
xsol=getvalue(x) # column T*1
fwsol=getvalue(fw)
aiwtsol=getvalue(aiwt)
return xsol,fwsol,aiwtsol
end
