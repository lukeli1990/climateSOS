using JuMP
using Gurobi
using MathProgBase
using MatpowerCases
using JLD
    # find out critical buses
    #run offline to determine the critical buses
    # Manually Changes: B jldname
     T=5
     B=14


     critbus=zeros(B,T)

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

     #modified 3 bus system by deleting load on bus 3
     loadval=mpc["bus"][:,3]*ones(1,T) #B*T

for busind=1:B
    for tind=1:T

 m2nd = Model(solver=GurobiSolver(LogToConsole=0))
 @variable(m2nd, Pgen[1:B] )
 @variable(m2nd, theta[1:B])
 @variable(m2nd, Pf[1:branchnum]) #14e power flow limit
 @variable(m2nd, aiwt[1:B],Bin)

 @constraint(m2nd,aiwt[busind]==0)

 #ref bus is the gen bus of lowest index; i.e. genbus[1]
 @constraint(m2nd, theta[1] == 0)

 #gen UB and LB 14(d)
 @constraint(m2nd, [i=1:B], Pgen[i]-aiwt[i]*Pmin[i]>=0)
 @constraint(m2nd, [i=1:B], -Pgen[i]+aiwt[i]*Pmax[i]>=0)

 #linelimit
 @constraint(m2nd,[i=1:branchnum], Pf[i]-linemin[i]*aiwt[fromto[i,1]]>=0)
 @constraint(m2nd,[i=1:branchnum], Pf[i]-linemin[i]*aiwt[fromto[i,2]]>=0)
 @constraint(m2nd,[i=1:branchnum],-Pf[i]+linemax[i]*aiwt[fromto[i,1]]>=0)
 @constraint(m2nd,[i=1:branchnum],-Pf[i]+linemax[i]*aiwt[fromto[i,2]]>=0)

 #14 (f) power balance ASSUME: load not changing for T (will change later)
 @constraint(m2nd, [i=1:B], edgeto[:,i]'*Pf[:]-edgefrom[:,i]'*Pf[:]+Pgen[i]-loadval[i,tind]==0)

 # 14 (g) is needed if any missing node is connecting with two buses
 @constraint(m2nd, [i=1:branchnum], -Pf[i] + (theta[fromto[i,1]]-theta[fromto[i,2]])/reacx[i]+BigM*(2-aiwt[fromto[i,1]]-aiwt[fromto[i,1]])>=0)
 @constraint(m2nd, [i=1:branchnum], Pf[i] - (theta[fromto[i,1]]-theta[fromto[i,2]])/reacx[i]-BigM*(aiwt[fromto[i,1]]+aiwt[fromto[i,1]]-2)>=0)



 status = solve(m2nd)

 if status == :Infeasible
 critbus[busind,tind]=1
 end

println("B=",busind,";T=",tind," is done")
end
end

save("criticalbus14.jld", "critbus", critbus)

## evaluate critical bus to see if they can support the system. If so this is the only combination
truecrit=zeros(1,T)
for t=1:T
 c2nd = Model(solver=GurobiSolver(LogToConsole=0))
 @variable(c2nd, Pgen[1:B] )
 @variable(c2nd, theta[1:B])
 @variable(c2nd, Pf[1:branchnum]) #14e power flow limit
 @variable(c2nd, aiwt[1:B],Bin)

 @constraint(c2nd,[i=1:B],aiwt[i]==critbus[i,t])

 #ref bus is the gen bus of lowest index; i.e. genbus[1]
 @constraint(c2nd, theta[1] == 0)

 #gen UB and LB 14(d)
 @constraint(c2nd, [i=1:B], Pgen[i]-aiwt[i]*Pmin[i]>=0)
 @constraint(c2nd, [i=1:B], -Pgen[i]+aiwt[i]*Pmax[i]>=0)

 #linelimit
 @constraint(c2nd,[i=1:branchnum], Pf[i]-linemin[i]*aiwt[fromto[i,1]]>=0)
 @constraint(c2nd,[i=1:branchnum], Pf[i]-linemin[i]*aiwt[fromto[i,2]]>=0)
 @constraint(c2nd,[i=1:branchnum],-Pf[i]+linemax[i]*aiwt[fromto[i,1]]>=0)
 @constraint(c2nd,[i=1:branchnum],-Pf[i]+linemax[i]*aiwt[fromto[i,2]]>=0)

 #14 (f) power balance ASSUME: load not changing for T (will change later)
 @constraint(c2nd, [i=1:B], edgeto[:,i]'*Pf[:]-edgefrom[:,i]'*Pf[:]+Pgen[i]-loadval[i,t]==0)

 # 14 (g) is needed if any missing node is connecting with two buses
 @constraint(c2nd, [i=1:branchnum], -Pf[i] + (theta[fromto[i,1]]-theta[fromto[i,2]])/reacx[i]+BigM*(2-aiwt[fromto[i,1]]-aiwt[fromto[i,1]])>=0)
 @constraint(c2nd, [i=1:branchnum], Pf[i] - (theta[fromto[i,1]]-theta[fromto[i,2]])/reacx[i]-BigM*(aiwt[fromto[i,1]]+aiwt[fromto[i,1]]-2)>=0)



 status2 = solve(c2nd)
 if status2 != :Infeasible
 truecrit[t]=1
 end

end

if minimum(truecrit)>0
println("This critbus is true in the sense that they alone can support the power system")
end
