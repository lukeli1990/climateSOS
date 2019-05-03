println()
println()
println()
println()
println()
println()
println("====================program start=========================")
using JuMP
using Gurobi
using Plots
include("first_stage.jl")
include("relax.jl")
include("upperbound.jl")
include("subproblem.jl")
include("second_stage.jl")
include("ranking.jl")

# ==========parameter initialization==================#
iter=1
itermax=2000


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




#============== order index of scenarios =======================#
orderscenw=zeros(Int,B,T,sizew)
orderscenw=rankw(scenariow)



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
count_fw_vio=zeros(itermax) #count the number of difference between 1st and 2nd fw
objv2nd=[0;sizew*ones(itermax-1)]


aiwt=zeros(B,T,sizew,itermax) #aiwt for all interations
index=1

#==========================main program========================#
while iter <= itermax && objv2nd[index]<3.9 # 20% 60%; 18.5 #90%






    #stage 1
	println("==============Stage 1: Iter:",iter,"==== START ================")


        (opttemp,fwtemp,aiwttemp) =f1st(x0,effcost,UB,scenariow,loadbus,loadnum,ei,di,iter,aiwt[:,:,:,1:index],orderscenw)
        stage1st_fw[:,iter]=fwtemp
		optsol[:,:,iter] = opttemp #B*T
        aiwt[:,:,:,iter] = aiwttemp #B*T*sizew
		objv1st[iter]=sum(opttemp[i,:]'*effcost[i,:] for i=1:B)


    println()
    println()
	println("==============Stage 2: Iter:",iter,"==== START ================")

	#======calculate the value of real fw  given x and a solved in the first stage  ==========================#
      fwtempreal=zeros(sizew)
      for scenind=1:sizew
	     fwtempreal[scenind]=subprob(aiwttemp[:,:,scenind],scenind)
	  end
      stage2nd_fw[:,iter]=fwtempreal
	  objv2nd[iter]=sum(fwtempreal)




    println()
    println()
	println("==============Stage 3: Iter:",iter,"==== START ================")

    count_fw=0
    for scenind=1:sizew
        if (fwtemp[scenind]==1 && fwtempreal[scenind]==0)
            count_fw=count_fw+1
        end
    end
    count_fw_vio[iter]=count_fw

    println()
    println()
    println("========================= All Stage Summary ===========================")
    println("==============Stage 1: Obj:",objv1st[iter],"; optimal solution:",optsol[:,:,iter],"===========")
    println("==============aiwt sum =",sum(aiwttemp)," ========= fw_1st sum=",sum(fwtemp))
    println("==============Stage 2: fw_2nd sum=",objv2nd[iter]," Iter: ",iter,"====================")


    iter=iter+1
	index=iter-1
end
