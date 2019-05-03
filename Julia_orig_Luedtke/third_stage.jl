function f3rd(a,Tk,Wk,bk)
    # No decomposition yet
    # a: B*T


    T=3
    #expand the sol of a
    anew=[a[:,1];a[:,2];a[:,3]] # [T=1;T=2;T=3]
    #14 bus param
    B=14 #bus number


    Tagg=cat([1,2],Tk,Tk,Tk) #  B*T cheng B*T
    Wagg=cat([1,2],Wk,Wk,Wk) #B*T cheng B*T
    bagg=[bk;bk;bk] #B*T


#dimension of pi
dimpi= size(Tagg,1)
sizey=size(Wk,2)

#=======feasibility check for Tk Wk bk=====#
# m3rd_fea = Model(solver=GurobiSolver(LogToConsole=0))
# @variable(m3rd_fea, y[1:sizey*T])
#
# @constraint(m3rd_fea,Tagg*anew+Wagg*y.>=bagg)
#
# status0 = solve(m3rd_fea)

#println("$status0")









#===========generate alpha and beta =============#


m3rd = Model(solver=GurobiSolver(LogToConsole=0))
@variable(m3rd, pidual[1:dimpi]>=0)

#ref bus is the gen bus of lowest index; i.e. genbus[1]
@constraint(m3rd, sum(pidual)==1)

#14 (f) power balance
@constraint(m3rd, pidual'*Wagg.==0)

@objective(m3rd, Max, sum(pidual[i]*(bagg[i]-Tagg[i,:]'*anew) for i=1:dimpi))

status1 = solve(m3rd)

piopt= getvalue(pidual)

alpha_cut= piopt'*Tagg # row vector

beta_cut= piopt'*bagg # scalar

#println("$status1")







#==================check that cut did separate solution from the feasible set==#
mfeas = Model(solver=GurobiSolver(LogToConsole=0))
 @variable(mfeas, afeas[1:B*T]) # a[T=1;T=2;T=3]
 @variable(mfeas, rec[1:sizey*T])# Pg[T=1;T=2;T=3]

 @objective(mfeas, Max,beta_cut - sum(alpha_cut[i]*afeas[i] for i=1:(B*T)))

 #ref bus is the gen bus of lowest index; i.e. genbus[1]
 @constraint(mfeas, Tagg*afeas+Wagg*rec.>=bagg)

# @constraint(mfeas,alpha_cut*afeas<=beta_cut) #orginally strict inequality

status2 = solve(mfeas)
 feasobjcost=getobjectivevalue(mfeas)
#
  if feasobjcost>=1e-5
  println("====================== ERROR: on generating alpha and beta==================")
  quit()
  end

return alpha_cut,beta_cut

end
