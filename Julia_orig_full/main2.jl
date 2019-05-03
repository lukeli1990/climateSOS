println()
println()
println()
println()
println()
println()
println("====================program start=========================")
using JuMP
using Gurobi
using MathProgBase

include("relax.jl")
include("upperbound.jl")


# ==========parameter initialization==================#
T=3 #time interval
B=14 #bus number
installcost=[1 4 3;2 1 5; 7 1 3;1 4 3;2 1 5; 7 1 3;1 4 3;2 1 5; 7 1 3;1 4 3;2 1 5; 7 1 3; 1 4 3; 2 1 5] # Bus* Time
effcost=[-3 1 3; 1 -4 5; 6 -2 3;-3 1 3; 1 -4 5; 6 -2 3;-3 1 3; 1 -4 5; 6 -2 3;-3 1 3; 1 -4 5; 6 -2 3;-3 1 3; 1 -4 5] # c1-c2,c2-c3,c3;B*T


# ==========14 bus system parameters==================#
gennum=5
genbus=[1; 2; 3; 6; 8]
Pmax=[362; 63 ; 20;0;0; 20;0; 20;0;0;0;0;0;0]
Pmin=[0 ;0; 0; 0; 0 ; 0;0;0;0;0;0;0;0;0]



branchnum=20
fromto=[1 2;1 5;2 3;2 4;2 5;3 4;4 5;4 7;4 9;5 6;6 11;6 12;6 13;7 8;7 9;9 10;9 14;10 11;12 13;13 14] #branchnum*2
edgeto=zeros(branchnum,B)#branchnum*busnum
edgefrom=zeros(branchnum,B)#branchnum*busnum
for i=1:branchnum
    edgeto[i,fromto[i,2]]=1
    edgeto[i,fromto[i,1]]=1
end
reacx=[0.06;0.2;0.2;0.2;0.2;0.2;0.04;0.2;0.6;0.3;0.1;0.3;0.2;0.2;0.1;0.1;0.3;0.2;0.2;0.4] #branch num
fmax=[480;130;150;160;160;160;670;150;60;120;140;110;200;170;270;330;100;140;100;80]
linemax=fmax
linemin=-fmax
BigM=sum((fmax[i]/reacx[i])   for i =1:branchnum)

#modified 3 bus system by deleting load on bus 3
loadbus=[2;3;4;5;6;9;10;11;12;13;14]
loadnum=11
loadval=[0;21;94;47;7;11;0;0;29;9;3;6;13;14] #busnum

#===============uncertainty initialization (random)===============#
sizew=20  #number of scenarios
srand(1234)  #fix seed for replication
scenariow=1+5*rand(B,T,sizew)
scenariow[:,1,:]=zeros(B,sizew)-1
ei=zeros(B)
di=ones(B)+[0;1;2;3;2;1;0;0;1;2;2;1;2;2]



#============== initial adaptation and upperbound=======================#
x0=[1;1;1;1;1;1;1;1;1;1;1;1;1;1]
UB=zeros(B)

for busind=1:B
    UB[busind]=upperboundx(scenariow[busind,:,:],ei[busind],di[busind])
end


#==========================save solutions========================#
optsol=zeros(B,T) #recorded for all iterations
optobjv=zeros(1)
optfw=zeros(sizew) #fw solved from first stage
optaiwt=zeros(B,T,sizew) #aiwt for all interations




#==========================main program========================#
    println()
    println()
	println("============== Optimization START ================")

    m1st = Model(solver=GurobiSolver(OutputFlag=0))

    #==========================first stage========================#
    @variable(m1st, x0[i]<=x[i=1:B,1:T]<= UB[i],Int)
    @variable(m1st, fw[1:sizew],Bin)
    @variable(m1st, aiwt[1:B,1:T,1:sizew],Bin)
    @variable(m1st, 0<=aiwtxit[i=1:B,1:T,1:sizew]<=UB[i]) #slack variable for product of awit and xit


    @objective(m1st, Min, sum(effcost[i,j]*x[i,j] for i=1:B,j=1:T))

    #initial condition = x0
    @constraint(m1st, [i=1:B], x[i,1]==x0[i])

    #increasing adaptation at each time step
    @constraint(m1st, [i=1:B,t=2:T], x[i,t]>=x[i,t-1])

    #sum of total fw
    @constraint(m1st, sum(fw[i] for i=1:sizew)>=8) #90 percent requirement 60 20 40


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



    #========== second stage constraints =================#

    @variable(m1st, Pgen[1:B,1:T,1:sizew])
    @variable(m1st, theta[1:B,1:T,1:sizew])
    @variable(m1st, Pf[1:branchnum,1:T,1:sizew]) #14e power flow limit

    #ref bus is the gen bus of lowest index; i.e. genbus[1]
    @constraint(m1st, [j=1:T,w=1:sizew],theta[1,j,w] == 0)

    #gen UB and LB 14(d)
    @constraint(m1st, [i=1:B,j=1:T,w=1:sizew], Pgen[i,j,w]-aiwt[i,j,w]*Pmin[i]>=0)
    @constraint(m1st, [i=1:B,j=1:T,w=1:sizew], -Pgen[i,j,w]+aiwt[i,j,w]*Pmax[i]>=0)

    #linelimit
    @constraint(m1st,[i=1:branchnum,j=1:T,w=1:sizew], Pf[i,j,w]-linemin[i]*aiwt[fromto[i,1],j,w]>=0)
    @constraint(m1st,[i=1:branchnum,j=1:T,w=1:sizew], Pf[i,j,w]-linemin[i]*aiwt[fromto[i,2],j,w]>=0)
    @constraint(m1st,[i=1:branchnum,j=1:T,w=1:sizew],-Pf[i,j,w]+linemax[i]*aiwt[fromto[i,1],j,w]>=0)
    @constraint(m1st,[i=1:branchnum,j=1:T,w=1:sizew],-Pf[i,j,w]+linemax[i]*aiwt[fromto[i,2],j,w]>=0)

    #14 (f) power balance ASSUME: load not changing for T (will change later)
    @constraint(m1st, [i=1:B,j=1:T,w=1:sizew], edgeto[:,i]'*Pf[:,j,w]-edgefrom[:,i]'*Pf[:,j,w]+Pgen[i,j,w]-loadval[i]>=(fw[w]-1)*BigM)
    @constraint(m1st, [i=1:B,j=1:T,w=1:sizew], edgeto[:,i]'*Pf[:,j,w]-edgefrom[:,i]'*Pf[:,j,w]+Pgen[i,j,w]-loadval[i]<=(1-fw[w])*BigM)

    # 14 (g) is needed if any missing node is connecting with two buses
    @constraint(m1st, [i=1:branchnum,j=1:T,w=1:sizew], -Pf[i,j,w] + (theta[fromto[i,1],j,w]-theta[fromto[i,2],j,w])/reacx[i]+BigM*(3-aiwt[fromto[i,1],j,w]-aiwt[fromto[i,1],j,w]-fw[w])>=0)
    @constraint(m1st, [i=1:branchnum,j=1:T,w=1:sizew], Pf[i,j,w] - (theta[fromto[i,1],j,w]-theta[fromto[i,2],j,w])/reacx[i]-BigM*(fw[w]+aiwt[fromto[i,1],j,w]+aiwt[fromto[i,1],j,w]-3)>=0)

    status=solve(m1st)
    println("Status = $status")

    optsol=getvalue(x)
    optobjv=getobjectivevalue(m1st)
    optaiwt=getvalue(aiwt)
    optfw=getvalue(fw)





	println("============== Objective cost:",optobjv,"; optimal solution:",optsol[:,:],"===========")
    println("==============aiwt sum =",sum(optaiwt),"; fw values:",optfw[:]," ========== ")
