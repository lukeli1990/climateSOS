function coeff_TWb()
#======== generating A matrix with 2nd stage problem============#
#======== Tx + Wy >= b ===========#
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
loadval=[0;21;94;47;7;11;0;0;29;9;3;6;13;14] #busnum

m2nd = Model(solver=GurobiSolver(LogToConsole=0))
@variable(m2nd, Pgen[1:B])
@variable(m2nd, theta[1:B])
@variable(m2nd, Pf[1:branchnum]) #14e power flow limit
@variable(m2nd, a[1:B])

#ref bus is the gen bus of lowest index; i.e. genbus[1]
@constraint(m2nd, theta[1] >= 0)
@constraint(m2nd, -theta[1]>= 0)

#gen UB and LB 14(d)
@constraint(m2nd, [i=1:B], Pgen[i]-a[i]*Pmin[i]>=0)
@constraint(m2nd, [i=1:B], -Pgen[i]+a[i]*Pmax[i]>=0)

#linelimit
@constraint(m2nd,[i=1:branchnum], Pf[i]-linemin[i]*a[fromto[i,1]]>=0)
@constraint(m2nd,[i=1:branchnum], Pf[i]-linemin[i]*a[fromto[i,2]]>=0)
@constraint(m2nd,[i=1:branchnum],-Pf[i]+linemax[i]*a[fromto[i,1]]>=0)
@constraint(m2nd,[i=1:branchnum],-Pf[i]+linemax[i]*a[fromto[i,2]]>=0)

#14 (f) power balance
@constraint(m2nd, [i=1:B], edgeto[:,i]'*Pf-edgefrom[:,i]'*Pf+Pgen[i]-loadval[i]>=0)
@constraint(m2nd, [i=1:B], -(edgeto[:,i]'*Pf-edgefrom[:,i]'*Pf+Pgen[i]-loadval[i])>=0)

# 14 (g) is needed if any missing node is connecting with two buses
@constraint(m2nd, [i=1:branchnum], -Pf[i] + (theta[fromto[i,1]]-theta[fromto[i,2]])/reacx[i]+BigM*(2-a[fromto[i,1]]-a[fromto[i,1]])>=0)
@constraint(m2nd, [i=1:branchnum], Pf[i] - (theta[fromto[i,1]]-theta[fromto[i,2]])/reacx[i]-BigM*(a[fromto[i,1]]+a[fromto[i,1]]-2)>=0)

JuMP.build(m2nd)
Amatr=MathProgBase.getconstrmatrix(internalmodel(m2nd))
linUB=MathProgBase.getconstrUB(internalmodel(m2nd))
linLB=MathProgBase.getconstrLB(internalmodel(m2nd))

for i=1:size(Amatr,1)
    if(linUB[i]!=Inf)
      println("Error on generating coefficient")
    end
end

T=Amatr[:,(end-B+1):end]
W=Amatr[:,1:(end-B)]
b=linLB


return T,W,b

end
