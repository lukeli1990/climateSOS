#presolve order 2
using JuMP
using Gurobi
using PowerModels
using JLD
    # find out critical buses and save power model parameters
    # run offline to determine the critical buses
    # Manually Changes: B jldname
     T=5
     B=14


     critbus=zeros(B,T)

     path="C:\\Users\\libow\\Desktop\\LANL_folder\\nesta\\nesta\\opf\\typ\\nesta_case14_ieee.m"
     mpc = PowerModels.parse_file(path)

     gennum=length(mpc["gen"])
     genbus=zeros(gennum,1)
     Pgmax=zeros(gennum,1)
     for i=1:gennum
     genbus[i]=mpc["gen"][string(i)]["gen_bus"] #size gen
     Pgmax[i]=mpc["gen"][string(i)]["pmax"]*100 #size gen;pu to MW
     end
     genbus=round.(Int, genbus)
     Pmax=zeros(B,1)
     Pmax[genbus]=Pgmax #save
     Pmin=zeros(B,1) #save


     branchnum=length(mpc["branch"])#save
     edgeto=zeros(Int,branchnum,B)#save
     edgefrom=zeros(Int,branchnum,B)#save
     fromto=zeros(Int,branchnum,2)#save
     reacx=zeros(branchnum)#save
     fmax=zeros(branchnum)
     for i=1:branchnum
         frombus=mpc["branch"][string(i)]["f_bus"]
         frombus=round(Int,frombus)
         tobus=mpc["branch"][string(i)]["t_bus"]
         tobus=round(Int,tobus)
         fromto[i,:]=[frombus tobus]
         edgeto[i,fromto[i,2]]=1
         edgefrom[i,fromto[i,1]]=1
         reacx[i]=mpc["branch"][string(i)]["br_x"]
         fmax[i]=mpc["branch"][string(i)]["rate_a"]*100 # pu to MW
     end
     linemax=fmax#save
     linemin=-fmax#save
     BigM=sum(fmax[i]   for i =1:branchnum)#save

     loadnum=length(mpc["load"])
     loadbus=zeros(Int,loadnum,1)
     loadvalbase=zeros(B,1)
     for i=1:loadnum
          loadbus=mpc["load"][string(i)]["load_bus"]
          loadbus=round(Int,loadbus)
          loadvalbase[loadbus]=mpc["load"][string(i)]["pd"]*100 #pu to MW
     end
     loadval=loadvalbase*ones(1,T) #save

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

println()
println()
if minimum(truecrit)>0
println("Crit bus support the system")
else
println("No critical bus support the system")
end

save("criticalbus14.jld", "critbus", critbus,"truecrit",truecrit)
save("powermodel14.jld", "Pmax", Pmax,"Pmin", Pmin,"branchnum",branchnum,"edgeto",edgeto,"edgefrom",edgefrom,"fromto",fromto,"reacx",
         reacx,"linemax",linemax,"linemin",linemin,"BigM",BigM,"loadval",loadval)
