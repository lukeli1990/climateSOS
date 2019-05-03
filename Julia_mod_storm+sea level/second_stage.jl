function f2nd(scenariow,aiwtscen)
# Manually Changes: B 

     T=5
     B=14


     mpc = loadcase("case14")

     gennum=size(mpc["gen"],1)
     genbus=round.(Int, mpc["gen"][:,1]) #size gen
     Pgmax=mpc["gen"][:,9]
     Pgmin=mpc["gen"][:,10]
     Pmax=zeros(B,1) #size B
     Pmax[genbus]=Pgmax
     Pmin=zeros(B,1) #size B
     Pmin[genbus]=Pgmin



     branchnum=size(mpc["branch"],1)
     fromto=round.(Int,mpc["branch"][:,1:2]) #integer
     edgeto=zeros(branchnum,B)#branchnum*busnum
     edgefrom=zeros(branchnum,B)#branchnum*busnum
     for i=1:branchnum
         edgeto[i,fromto[i,2]]=1
         edgeto[i,fromto[i,1]]=1
     end
     reacx=mpc["branch"][:,4] #branch num
    # fmax=mpc["branch"][:,6] #it is zero in data
     fmax=[480;130;150;160;160;160;670;150;60;120;140;110;200;170;270;330;100;140;100;80]
     linemax=fmax
     linemin=-fmax
     BigM=sum(fmax[i]   for i =1:branchnum)

     loadval=mpc["bus"][:,3]*ones(1,T) #B*T


 m2nd = Model(solver=GurobiSolver(LogToConsole=0))
 @variable(m2nd, Pgen[1:B,1:T])
 @variable(m2nd, theta[1:B,1:T])
 @variable(m2nd, Pf[1:branchnum,1:T])
 @variable(m2nd, aiwt[1:B,1:T]) #does not need to be binary


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
