# workspace()
# Settings: change into 14 bus system for double check
# Simulate full problem with original scenario
# Manually Changes: B sizew jldname
println()
println()
println()
println()
println()
println()
println("====================program start=========================")
using JuMP
using Gurobi
using MatpowerCases
using JLD


# ==========main parameter initialization==================#
T=5
sizew=2000
B=14
conflvl=0.1 #1-epsilon

# ==========hardening parameters========================#
#Note: current elevation are manually set
ei=ones(B,T) #Current Elevation of the Bus (worst case elevation)




# ===========cost parameter initialization=================#
# Note: current cost are randomly generated
srand(1234)  #fix seed for replication
installcost=rand(1:10,B,T) #B*T
effcost=[installcost[:,1]-installcost[:,2] installcost[:,2]-installcost[:,3] installcost[:,3]-installcost[:,4] installcost[:,4]-installcost[:,5] installcost[:,5]] #c1-c2,c2-c3,c3;B*T



#===============scenario initialization===============#
#readin scenario
temp2=load("14scenario2000.jld")
#temp2=load("scenario20.jld")
scenarioraw=temp2["scenarior"] #B * T * sizew
println()
println()
println("====================Data Reading Finish=========================")

#modified scenario: scenario-e
scenariow=zeros(B,T,sizew) #B * T * sizew
for i=1:B
    for t=1:T
        for s=1:sizew
            scenariow[i,t,s]=scenarioraw[i,t,s]-ei[i,t]
        end
    end
end



#generate UB for h
trivialUB=zeros(B,T)
for i=1:B
    for t=1:T
        trivialUB[i,t]=max(maximum(scenariow[i,t,:]),0)
    end
end
for i=1:B
    for t=2:T
        trivialUB[i,t]=max(trivialUB[i,t],trivialUB[i,t-1])
    end
end




#==========================save solutions========================#
hsol=zeros(B,T) #recorded for all iterations
fwsol=zeros(sizew) #fw solved from first stage
aiwtsol=zeros(B,T,sizew) #aiwt for all iterations: scenariosort



#==========================power parameters solutions========================#
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




tic()
#==========================main program========================#
m1st = Model(solver=GurobiSolver())


@variable(m1st, 0<=h[i=1:B,t=1:T]<=trivialUB[i,t])
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
@constraint(m1st, [i=1:B,t=1:T,s=1:sizew], maximum(scenariow)*(1-fw[s]) + h[i,t] >= scenariow[i,t,s]*aiwt[i,t,s]) #90 percent requirement 60 20 40



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
