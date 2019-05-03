#presolve order 1
#run offline read scenario from full scenario sets
# Manually Changes: B sizew jldname folderid
using JSON
using JLD

B=14
T=5
sizew=100


scenarior=zeros(B,T,sizew)
SSscen=zeros(B,T,sizew)
SLscen=zeros(B,T,sizew)

s = open("C:\\Users\\libow\\Desktop\\LANL_folder\\Julia_SOS1\\with sea lvl and ld redistributin\\SOS1_14_single_run\\nf2_all.json") do file
    read(file, String)
end

data=JSON.parse(s)


indrand=rand(1:5000,1,sizew)
indbus=rand(1:118,1,B)

for s1=1:sizew
    for t=1:T
           SLscen[:,t,s1]=data["SL"][string(indrand[s1])][t]*ones(B,1)

        for i=1:B
           SSscen[i,t,s1]=data["SS"][string(indrand[s1])][t][indbus[i]]
        end
    end
end

eitot = [1,1,1,6,4,1,4,8,4,1,6,8,8,6,6,4,5,6,7,5,4,1,1,4,2,4,2,1,1,6,2,3,6,3,1,2,2,3,1,2,2,2,1,1,2,2,2,3,3,4,5,2,
        1,
        5,
        1,
        2,
        3,
        2,
        1,
        1,
        2,
        3,
        1,
        2,
        2,
        3,2,3,2,2,4,3,5,2,2,1,2,1,1,3,1,1,2,2,4,2,5,7,4,3,5,5,5,1,3,2,1,2,5,2,2,3,3,2,4,1,3,4,5,5,6,6,4,1,1,3,1,1]# for 118 bus
        ei=eitot[indbus]#1*B

        scenarior=SLscen+SSscen

save("14scenario100+sea.jld", "scenarior", scenarior,"SLscen",SLscen,"ei",ei)

if minimum(ei) - maximum(SLscen)<0
   println()
   println()
   println("Next step: need to run load redistribution")
end
