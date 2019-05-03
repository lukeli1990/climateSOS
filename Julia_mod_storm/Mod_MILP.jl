# workspace()
# Settings: change into 14 bus system for double check
# Simulate full problem with modified scenario
# Manually Changes: B sizew jldname Int/not
tic()
println()
println()
println()
println()
println()
println()
println("====================program start=========================")
using JuMP
using Gurobi
using PowerModels
using JLD


# ==========main parameter initialization==================#
T=5
sizew=2000
B=14
conflvl=0.1 #1-epsilon

# ==========hardening parameters========================#
#readin scenario
temp2=load("14scenario2000.jld")
ei=temp2["ei"]'*ones(1,T) #Current Elevation of the Bus




# ===========cost parameter initialization=================#
# Note: current cost are randomly generated
srand(1234)  #fix seed for replication
installcost=rand(800:1200,B,T) #B*T
effcost=[installcost[:,1]-installcost[:,2] installcost[:,2]-installcost[:,3] installcost[:,3]-installcost[:,4] installcost[:,4]-installcost[:,5] installcost[:,5]] #c1-c2,c2-c3,c3;B*T



#===============scenario initialization===============#
scenarioraw=temp2["scenarior"] #B * T * sizew
println()
println()
println("====================Data Reading Finish=========================")

#modified scenario: max(scenario-e,0)
scenariow=zeros(B,T,sizew) #B * T * sizew
for i=1:B
    for t=1:T
        for s=1:sizew
            scenariow[i,t,s]=max(scenarioraw[i,t,s]-ei[i,t],0)
        end
    end
end



#generate UB for h
trivialUB=zeros(B)
for i=1:B
trivialUB[i]=maximum(scenariow[i,:,:])
end
trivialUB=ceil.(trivialUB)#Int/not



#==========================save solutions========================#
hsol=zeros(B,T) #recorded for all iterations
fwsol=zeros(sizew) #fw solved from first stage
aiwtsol=zeros(B,T,sizew) #aiwt for all iterations: scenariosort



#==========================power parameters solutions========================#
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





#==========================main program========================#
m1st = Model(solver=GurobiSolver())


@variable(m1st, 0<=h[i=1:B,t=1:T]<=trivialUB[i],Int)#Int/not
@variable(m1st, fw[1:sizew],Bin)#sizew is indexed with orignal order
@variable(m1st, aiwt[1:B,1:T,1:sizew], Bin) #sizew is indexed with scenario ascending order
@variable(m1st, Pgen[1:B,1:T,1:sizew])
@variable(m1st, theta[1:B,1:T,1:sizew])
@variable(m1st, Pf[1:branchnum,1:T,1:sizew])





@objective(m1st, Min, sum(effcost[i,t]*h[i,t] for i=1:B,t=1:T))


#increasing adaptation at each time step
@constraint(m1st, [i=1:B,t=2:T], h[i,t]>=h[i,t-1])



#sum of total fw
@constraint(m1st, sum(fw[s] for s=1:sizew)>=sizew*conflvl) #90 percent requirement 60 20 40


#second stage
@constraint(m1st, [i=1:B,t=1:T,s=1:sizew], scenariow[i,t,s]*(1-fw[s]) + h[i,t] >= scenariow[i,t,s]*aiwt[i,t,s]) #90 percent requirement 60 20 40



 #gen UB and LB 14(d)
 @constraint(m1st, [i=1:B,t=1:T,s=1:sizew], Pgen[i,t,s]-aiwt[i,t,s]*Pmin[i]>=0)
 @constraint(m1st, [i=1:B,t=1:T,s=1:sizew], -Pgen[i,t,s]+aiwt[i,t,s]*Pmax[i]>=0)

 #linelimit
 @constraint(m1st,[i=1:branchnum,t=1:T,s=1:sizew], Pf[i,t,s]-linemin[i]*aiwt[fromto[i,1],t,s]>=0)
 @constraint(m1st,[i=1:branchnum,t=1:T,s=1:sizew], Pf[i,t,s]-linemin[i]*aiwt[fromto[i,2],t,s]>=0)
 @constraint(m1st,[i=1:branchnum,t=1:T,s=1:sizew],-Pf[i,t,s]+linemax[i]*aiwt[fromto[i,1],t,s]>=0)
 @constraint(m1st,[i=1:branchnum,t=1:T,s=1:sizew],-Pf[i,t,s]+linemax[i]*aiwt[fromto[i,2],t,s]>=0)

 #14 (f) power balance ASSUME: load not changing for T (will change later)
 @constraint(m1st, [i=1:B,t=1:T,s=1:sizew], edgeto[:,i]'*Pf[:,t,s]-edgefrom[:,i]'*Pf[:,t,s]+Pgen[i,t,s]-loadval[i,t]==0)

 # 14 (g) is needed if any missing node is connecting with two buses
 @constraint(m1st, [i=1:branchnum,t=1:T,s=1:sizew], -Pf[i] + (theta[fromto[i,1],t,s]-theta[fromto[i,2],t,s])/reacx[i]+BigM*(2-aiwt[fromto[i,1],t,s]-aiwt[fromto[i,1],t,s])>=0)
 @constraint(m1st, [i=1:branchnum,t=1:T,s=1:sizew], Pf[i] - (theta[fromto[i,1],t,s]-theta[fromto[i,2],t,s])/reacx[i]-BigM*(aiwt[fromto[i,1],t,s]+aiwt[fromto[i,1],t,s]-2)>=0)



status = solve(m1st)
hsol=getvalue(h)
fwsol=getvalue(fw)
aiwtsol=getvalue(aiwt)
objcost=getobjectivevalue(m1st)
println("====================program finish successfully with Obj:",objcost,"; sum of fw=",sum(fwsol),"=====================")

toc()
