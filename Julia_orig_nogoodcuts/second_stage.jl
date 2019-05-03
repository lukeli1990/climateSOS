function f2nd(a) # a in COLUMN
    # decompose over scenario then time
    # a: busnum*1 column

    #3 bus param
    B=14 #bus number

    gennum=5
    genbus=[1; 2; 3; 6; 8]
    Pmax=[362; 63 ; 20; 20; 20].*a[genbus]
    Pmin=[0 ;0; 0; 0; 0].*a[genbus]
    genmat=zeros(14,5)
    genmat[1,1]=1
    genmat[2,2]=1
    genmat[3,3]=1
    genmat[6,4]=1
    genmat[8,5]=1


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
    linemax=a[fromto[:,1]].*a[fromto[:,2]].*fmax
    linemin=-linemax

    BigM=sum((fmax[i]/reacx[i])   for i =1:branchnum)
    BigMmin=BigM*(a[fromto[:,1]].*a[fromto[:,2]]-1)
    BigMmax=-BigMmin

    #modified 3 bus system by deleting load on bus 3
    loadbus=[2;3;4;5;6;9;10;11;12;13;14]
    loadnum=11
    loadval=[0;21;94;47;7;11;0;0;29;9;3;6;13;14] #busnum

m2nd = Model(solver=GurobiSolver(LogToConsole=0))
@variable(m2nd, Pmin[i] <= Pgen[i=1:gennum] <= Pmax[i])
@variable(m2nd, theta[1:B])
@variable(m2nd, linemin[i]<=Pf[i=1:branchnum]<=linemax[i]) #14e power flow limit

@expression(m2nd, Pfeff[i=1:branchnum], Pf[i] - (theta[fromto[i,1]]-theta[fromto[i,2]])/reacx[i])

#ref bus is the gen bus of lowest index; i.e. genbus[1]
@constraint(m2nd, theta[1] == 0)

#14 (f) power balance
@constraint(m2nd, [i=1:B], edgeto[:,i]'*Pf-edgefrom[:,i]'*Pf+genmat[i,:]'*Pgen-loadval[i]==0)

# 14 (g) is needed if any missing node is connecting with two buses
@constraint(m2nd, [i=1:branchnum], Pfeff[i]<=BigMmax[i])
@constraint(m2nd, [i=1:branchnum], BigMmin[i]<=Pfeff[i])
status = solve(m2nd)

println("$status")

if status == :Optimal
return 1
else
return 0
end
end
