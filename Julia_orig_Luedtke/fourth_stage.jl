function f4th(alpha,x0,UB,scenariow) #
    # No decomposition yet
    # alpha: 1 * (B*T); x0 B; UB B; scenariow B*T


    T=3
    #14 bus param
    B=14 #bus number

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
    loadval=[0 ;21; 94; 47; 7 ;11; 0; 0; 29; 9; 3; 6; 13; 14]*ones(1,T) #B*T


#===========generate alpha and beta =============#


m4th = Model(solver=GurobiSolver(OutputFlag=0))


@variable(m4th, x0[i]<=x[i=1:B,1:T]<= UB[i],Int)
@variable(m4th, aiwt[1:B,1:T],Bin)
@variable(m4th, 0<=aiwtxit[i=1:B,1:T]<=UB[i]) #slack variable for product of awit and xit

@objective(m4th, Min, sum(alpha[i+B*(j-1)]*aiwt[i,j] for i=1:B,j=1:T))

#initial condition = x0
@constraint(m4th, [i=1:B], x[i,1]==x0[i])

#increasing adaptation at each time step
@constraint(m4th, [i=1:B,t=2:T], x[i,t]>=x[i,t-1])


#constraint 13a with slack
@constraint(m4th, [i=1:B,j=1:T], (2*aiwt[i,j]-1)*scenariow[i,j] + ei[i]+di[i]*x[i,j] - 2*di[i]*aiwtxit[i,j] <=0)


#mccormick for aiwtxit
for i=1:B
    for j=1:T

            mccormick(m4th,aiwtxit[i,j],x[i,j],aiwt[i,j],x0[i],UB[i],0,1)

    end
end


    @constraint(m4th, [i=1:loadnum,j=1:T], aiwt[loadbus[i],j]>=1)



#========== second stage constraints =================#

@variable(m4th, Pgen[1:B,1:T])
@variable(m4th, theta[1:B,1:T])
@variable(m4th, Pf[1:branchnum,1:T]) #14e power flow limit

#ref bus is the gen bus of lowest index; i.e. genbus[1]
@constraint(m4th, [j=1:T],theta[1,j] == 0)

#gen UB and LB 14(d)
@constraint(m4th, [i=1:B,j=1:T], Pgen[i,j]-aiwt[i,j]*Pmin[i]>=0)
@constraint(m4th, [i=1:B,j=1:T], -Pgen[i,j]+aiwt[i,j]*Pmax[i]>=0)

#linelimit
@constraint(m4th,[i=1:branchnum,j=1:T], Pf[i,j]-linemin[i]*aiwt[fromto[i,1],j]>=0)
@constraint(m4th,[i=1:branchnum,j=1:T], Pf[i,j]-linemin[i]*aiwt[fromto[i,2],j]>=0)
@constraint(m4th,[i=1:branchnum,j=1:T],-Pf[i,j]+linemax[i]*aiwt[fromto[i,1],j]>=0)
@constraint(m4th,[i=1:branchnum,j=1:T],-Pf[i,j]+linemax[i]*aiwt[fromto[i,2],j]>=0)

#14 (f) power balance ASSUME: load not changing for T (will change later)
@constraint(m4th, [i=1:B,j=1:T], edgeto[:,i]'*Pf[:,j]-edgefrom[:,i]'*Pf[:,j]+Pgen[i,j]-loadval[i,j]==0)

# 14 (g) is needed if any missing node is connecting with two buses
@constraint(m4th, [i=1:branchnum,j=1:T], -Pf[i,j] + (theta[fromto[i,1],j]-theta[fromto[i,2],j])/reacx[i]+BigM*(2-aiwt[fromto[i,1],j]-aiwt[fromto[i,1],j])>=0)
@constraint(m4th, [i=1:branchnum,j=1:T], Pf[i,j] - (theta[fromto[i,1],j]-theta[fromto[i,2],j])/reacx[i]-BigM*(aiwt[fromto[i,1],j]+aiwt[fromto[i,1],j]-2)>=0)

status=solve(m4th)
if status != :Optimal
    println("====================== ERROR: on generating hvalue==================")
    quit()
end
#println("$status")
return getobjectivevalue(m4th),getvalue(aiwt)

end
