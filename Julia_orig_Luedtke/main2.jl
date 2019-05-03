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
include("first_stage.jl")
include("relax.jl")
include("upperbound.jl")
include("subproblem.jl")
include("second_stage.jl")
include("third_stage.jl")
include("fourth_stage.jl")
include("coeff_T_W_b_.jl")
include("LB_stage.jl")
#include("cut_generation_simple.jl")

# ==========parameter initialization==================#
iter=1
itermax=10


T=3 #time interval
B=14 #bus number

# modified 3 bus system bus 3 does not have load anymore
loadbus=[2;3;4;5;6;9;10;11;12;13;14]
loadnum=11


installcost=[1 4 3;2 1 5; 7 1 3;1 4 3;2 1 5; 7 1 3;1 4 3;2 1 5; 7 1 3;1 4 3;2 1 5; 7 1 3; 1 4 3; 2 1 5] # Bus* Time
effcost=[-3 1 3; 1 -4 5; 6 -2 3;-3 1 3; 1 -4 5; 6 -2 3;-3 1 3; 1 -4 5; 6 -2 3;-3 1 3; 1 -4 5; 6 -2 3;-3 1 3; 1 -4 5] # c1-c2,c2-c3,c3;B*T



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
optsol=zeros(B,T,itermax) #recorded for all iterations
objv1st=zeros(itermax)
stage1st_fw=zeros(sizew,itermax) #fw solved from first stage
stage2nd_fw=zeros(sizew,itermax) #fw solved from second stage
stage2nd_T=zeros(sizew,itermax) #Time index when fw fails
count_fw_vio=zeros(Int,itermax) #count the number of difference between 1st and 2nd fw
fw_vio_ind=zeros(Int,sizew,itermax) #record which scenario create difference between 1st and 2nd stage
objv2nd=[0;sizew*ones(Int,itermax-1)]

alpha_col=zeros(B*T,sizew,itermax)
beta_val=zeros(sizew,itermax)
hvalue=zeros(sizew,itermax) #vionum
LBvalue=zeros(sizew,itermax) #vionum

aiwt=zeros(B,T,sizew,itermax) #aiwt for all interations
index=1

conflvl=0.4 #90% or 60% or 20%  40%
pint=floor(Int,(1-conflvl)*sizew) #type integer


#====================== Coefficient matrix T W b===================#
(Tk,Wk,bk)=coeff_TWb() #the coeffcient will change if dit changes over time will change in future



#==========================main program========================#
while iter <= itermax && objv2nd[index]<ceil(Int,conflvl*sizew) #size20 90%=18 60%=12






    #stage 1
	println("==============Stage 1: Iter:",iter,"==== START ================")


        (opttemp,fwtemp,aiwttemp) =f1st(x0,effcost,UB,scenariow,loadbus,loadnum,ei,di,iter,count_fw_vio,alpha_col,fw_vio_ind,hvalue,LBvalue,conflvl)
        stage1st_fw[:,iter]=fwtemp
		optsol[:,:,iter] = opttemp #B*T
        aiwt[:,:,:,iter] = aiwttemp #B*T*sizew
		objv1st[iter]=sum(opttemp[i,:]'*effcost[i,:] for i=1:B)
        xvec=[opttemp[:,1];opttemp[:,2];opttemp[:,3]] #vectorize x for debug


	println("==============Stage 1: Obj:",objv1st[iter],"; optimal solution:",optsol[:,:,iter],"===========")
    println("==============aiwt sum =",sum(aiwttemp)," ========= fw_1st sum=",sum(fwtemp))


    println()
    println()
	println("==============Stage 2: Iter:",iter,"==== START ================")

	#======calculate the value of real fw  given x and a solved in the first stage  ==========================#
      fwtempreal=zeros(sizew)
      for scenind=1:sizew
	     (fwtempreal[scenind],stage2nd_T[scenind,iter])=subprob(aiwttemp[:,:,scenind],scenind)
	  end
      stage2nd_fw[:,iter]=fwtempreal
	  objv2nd[iter]=sum(fwtempreal)


	println("==============Stage 2: fw_2nd sum=",objv2nd[iter]," Iter: ",iter,"====================")

    println()
    println()
	println("==============Stage 3: Iter:",iter,"==== START ================")

    count_fw=0
    for scenind=1:sizew
        if (fwtemp[scenind]==1 && fwtempreal[scenind]==0) #identify violation index
            count_fw=count_fw+1
           (alpha_cut,beta_cut)=f3rd(aiwt[:,:,scenind,iter],Tk,Wk,bk)
           alpha_col[:,count_fw,iter]=alpha_cut
           beta_val[count_fw,iter]=beta_cut
           fw_vio_ind[count_fw,iter]=scenind
           #println(alpha_cut) #for debugging
        end
    end
    count_fw_vio[iter]=count_fw

    println()
    println()
    println("==============Stage 3: Iter:",iter,"======= Alpha value calculated =============")
           atest=zeros(B,T)# for debug
       for vioind=1:count_fw
           alpha_h=alpha_col[:,vioind,iter]
          (hvalue[vioind,iter],atest)=f4th(alpha_h,x0,UB,scenariow[:,:,fw_vio_ind[vioind,iter]])
          LBvalue[vioind,iter]=fLB(alpha_h)
           #println(checkval) #for debugging
       end
    #
      println()
      println()
      println("==============Stage 3: Iter:",iter,"======= Hcost and LB value calculated =============")
    #
    #

     iter=iter+1
	 index=iter-1
end
