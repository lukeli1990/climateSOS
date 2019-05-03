# workspace()
# Settings: change into 14 bus system for double check + evaluate Luedtke's cuts
# Methods: prf SOS1 + prf cont 2nd stage + prf/alg:Luedtke cut (opt/converg/stronger) + UB +
#          alg:critical buses from LB technique + prf:scenario dominance + modification on overall algorithm/discussion on previous contribution
# system: hardening + stormsurge +  load redistribution + load controllability + sea level rise + sum load threshold
# future work (new things): with expansion + new method (TSMILP or scenario reduction) + special warm start + h LB/UB + smarter constraints
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
using MathProgBase
using MatpowerCases
using JLD
include("scenario_dom.jl")
include("scenario_manipulate.jl")
include("first_stage.jl")
include("second_stage.jl")


# ==========main parameter initialization==================#
T=5
sizew=2000
B=14
conflvl=0.3 #1-epsilon

# ==========hardening parameters========================#
#Note: current elevation are manually set
ei=ones(B,T) #Current Elevation of the Bus (worst case elevation)
#Identify the critical buses that can't afford to lose
#Note: the current load values does not change with time
temp1=load("criticalbus14.jld")
critbus=temp1["critbus"]



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

#modified scenario: max(scenario-e,0)
scenariow=zeros(B,T,sizew) #B * T * sizew
for i=1:B
    for t=1:T
        for s=1:sizew
            scenariow[i,t,s]=max(scenarioraw[i,t,s]-ei[i,t],0)
        end
    end
end

#scenario dominance
domscen=scendom(scenariow) #sizew * sizew; for i!=j domscen(i,j)=1 ---> scen i >= scen j

#scenario manipulation
(thresind,rankscen,scenariosort)=scenmanipulate(scenariow, conflvl, critbus)
#remark: the rankscen contains index corresponding entries in scenariow; scenariosort[i]=scenariow[rankscen[i]]

#reverse index
invrankscen=zeros(Int,B,T,sizew)
for i=1:B
    for t=1:T
        invrankscen[i,t,:]=sortperm(rankscen[i,t,:])
    end
end
#remark: the invrankscen contains index corresponding entries in scenariosort; scenariow[i]=scenariosort[invrankscen[i]]

#generate UB for h
trivialUB=zeros(B,T)
trivialUB[:,1]=scenariosort[:,1,sizew]
for i=1:B
    for t=2:T
        trivialUB[i,t]=max(scenariosort[i,t,sizew],trivialUB[i,t-1])
    end
end

#index matrix track down which entry of threshold need SOS1
numSOS=0
SOSind=zeros(B,T)
SOSsize=zeros(B,T)
for i=1:B
    for t=1:T
        if thresind[i,t]<sizew
           numSOS=numSOS+1
           SOSind[i,t]=numSOS
           SOSsize[i,t]=sizew-thresind[i,t]+1 #include zero as first entry
        end
    end
end




#==========================save solutions========================#
hsol1st=zeros(B,T) #recorded for all iterations
objv1st=0# objective for all iterations
fw1st=zeros(1,sizew) #fw solved from first stage
fw2nd=zeros(1,sizew) #fw solved from second stage
fwsum2nd=0 #sum of fw from second stage
aiwt1st=zeros(B,T,sizew) #aiwt for all iterations: scenariosort
sos1st=zeros(B,T,sizew) #sos for all iterations: scenariosort
aiwtscen=zeros(B,T,sizew) #aiwt for all iterations: scenariow




tic()
#==========================main program========================#

println()
println()
println("==============Stage 1: START =====================")
       (hsol1st,sos1st,fw1st,aiwt1st,objv1st)=f1st(SOSsize,thresind,trivialUB,rankscen,scenariosort,conflvl,critbus,domscen,effcost)
println("==============Stage 1: Done with Obj:",objv1st,"===============")

#reformulate aiwt1st back to aiwtscen associated with scenariow
for i=1:B
    for t=1:T
        for s=1:sizew
        aiwtscen[i,t,s]=aiwt1st[i,t,invrankscen[i,t,s]]
        end
    end
end


println()
println()
println("==============Stage 2: START ========================")
#======calculate the value of real fw  given h and aUB solved in the first stage using scenario dominance  ==========================#
countcheck=0
domscenself=domscen+eye(sizew)
for s=1:sizew #s for scenariow
     if fw2nd[s]==0
     fw2nd[s]=f2nd(scenariow[:,:,s],aiwtscen[:,:,s])
     countcheck=countcheck+1
        if fw2nd[s]==1
            for s2=1:sizew
                fw2nd[s2]=max(fw2nd[s2],domscenself[s,s2]*fw2nd[s])
            end
        end
     end
end

fwsum2nd=sum(fw2nd)

println("=====Stage 2: Done with fw_sum=",fwsum2nd,"; larger than threshold: ",fwsum2nd>=conflvl*sizew,"; numofcheck:",countcheck,"/",sizew,"==============")

if fwsum2nd>=conflvl*sizew
    println()
    println()
    println("====================program finish successfully with Obj:",objv1st,"=====================")
end
toc()
