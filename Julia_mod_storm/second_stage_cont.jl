function f2nd(scenariow,aiwtscen,mpcmodel)
# Manually Changes: B

     T=5
     B=14


     Pmax=mpcmodel["Pmax"]
     Pmin=mpcmodel["Pmin"]
     branchnum=mpcmodel["branchnum"]
     edgeto=mpcmodel["edgeto"]
     edgefrom=mpcmodel["edgefrom"]
     fromto=mpcmodel["fromto"]
     reacx=mpcmodel["reacx"]
     linemax=mpcmodel["linemax"]
     linemin=mpcmodel["linemin"]
     BigM=mpcmodel["BigM"]
     loadval=mpcmodel["loadval"]


 m2nd = Model(solver=GurobiSolver(LogToConsole=0))
 @variable(m2nd, Pgen[1:B,1:T])
 @variable(m2nd, theta[1:B,1:T])
 @variable(m2nd, Pf[1:branchnum,1:T])
 @variable(m2nd, aiwt[1:B,1:T]) #does not need to be binary;for ModR_SOS


 #ref bus is the gen bus of lowest index; i.e. genbus[1]
 @constraint(m2nd, theta[1] == 0)

 #UB and LB for aiwt
 @constraint(m2nd, [i=1:B,t=1:T], aiwt[i,t]>=0)
 @constraint(m2nd, [i=1:B,t=1:T], -aiwt[i,t]>=-aiwtscen[i,t])

 #gen UB and LB 14(d)
 @constraint(m2nd, [i=1:B,t=1:T], Pgen[i,t]-aiwt[i,t]*Pmin[i]>=0)
 @constraint(m2nd, [i=1:B,t=1:T], -Pgen[i,t]+aiwt[i,t]*Pmax[i]>=0)

 #linelimit
 @constraint(m2nd,[i=1:branchnum,t=1:T], Pf[i,t]-linemin[i]*aiwt[fromto[i,1],t]>=0)
 @constraint(m2nd,[i=1:branchnum,t=1:T], Pf[i,t]-linemin[i]*aiwt[fromto[i,2],t]>=0)
 @constraint(m2nd,[i=1:branchnum,t=1:T],-Pf[i,t]+linemax[i]*aiwt[fromto[i,1],t]>=0)
 @constraint(m2nd,[i=1:branchnum,t=1:T],-Pf[i,t]+linemax[i]*aiwt[fromto[i,2],t]>=0)

 #14 (f) power balance ASSUME: load not changing for T (will change later)
 @constraint(m2nd, [i=1:B,t=1:T], edgeto[:,i]'*Pf[:,t]-edgefrom[:,i]'*Pf[:,t]+Pgen[i,t]-loadval[i,t]==0)

 # 14 (g) is needed if any missing node is connecting with two buses
 @constraint(m2nd, [i=1:branchnum,t=1:T], -Pf[i] + (theta[fromto[i,1],t]-theta[fromto[i,2],t])/reacx[i]+BigM*(2-aiwt[fromto[i,1],t]-aiwt[fromto[i,1],t])>=0)
 @constraint(m2nd, [i=1:branchnum,t=1:T], Pf[i] - (theta[fromto[i,1],t]-theta[fromto[i,2],t])/reacx[i]-BigM*(aiwt[fromto[i,1],t]+aiwt[fromto[i,1],t]-2)>=0)



 status = solve(m2nd)

 if status == :Infeasible
 return 0
 else
 return 1
 end

end
