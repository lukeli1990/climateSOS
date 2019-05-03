function f1st(SOSsize,thresind,trivialUB,rankscen,scenariosort,conflvl,critbus,domscen,effcost)
    # Manually Changes: B sizew Int/not

    T=5
    B=14
    sizew=2000

#m1st = Model(solver=GurobiSolver(OutputFlag=0))
m1st = Model(solver=GurobiSolver())


@variable(m1st, 0<=h[i=1:B,t=1:T]<=trivialUB[i],Int)#Int/not
@variable(m1st, 0<=sos[1:B,1:T,1:sizew]<=1) #sizew is indexed with scenario ascending order
@variable(m1st, fw[1:sizew],Bin)#sizew is indexed with orignal order
@variable(m1st, 0<=aiwtUB[1:B,1:T,1:sizew]<=1) #sizew is indexed with scenario ascending order

#define SOS and aiwt
for i=1:B
  for t=1:T
      if SOSsize[i,t]==0
          @constraint(m1st, [s=1:sizew], sos[i,t,s]==0)
          @constraint(m1st, [s=1:sizew], aiwtUB[i,t,s]==1)
      else
          #for sos[i,t,1:thresind[i,t]-1], we don't have requirement
          @constraint(m1st, [s=1:thresind[i,t]], aiwtUB[i,t,s]==1)
          addSOS1(m1st, sos[i,t,thresind[i,t]:sizew])
          @constraint(m1st, sum(sos[i,t,s] for s=thresind[i,t]:sizew)==1)
          @constraint(m1st, [s=(thresind[i,t]+1):sizew], aiwtUB[i,t,s]==sum(sos[i,t,s2] for s2=s:sizew))
      end
  end
end



@objective(m1st, Min, sum(effcost[i,t]*h[i,t] for i=1:B,t=1:T))


#increasing adaptation at each time step
@constraint(m1st, [i=1:B,t=2:T], h[i,t]>=h[i,t-1])


#relationship between h and sos
 for i=1:B
     for t=1:T
         @constraint(m1st, h[i,t]>=sum(scenariosort[i,t,s]*sos[i,t,s] for s=thresind[i,t]:sizew))
     end
 end


#sum of total fw
@constraint(m1st, sum(fw[s] for s=1:sizew)>=sizew*conflvl) #90 percent requirement 60 20 40


#critical buses
   for i=1:B
       for t=1:T
           if critbus[i,t]==1
              @constraint(m1st, [s=1:sizew], aiwtUB[i,t,s]>=fw[rankscen[i,t,s]])
          end
       end
   end

#scenario dominance
  for s1=1:sizew
     for s2=1:sizew
          if(domscen[s1,s2]==1)
              @constraint(m1st, fw[s1] <= fw[s2])
          end
      end
  end


#Luedtke cuts



status = solve(m1st)
println("$status")
if status == :Infeasible
println("============= First Stage Infeasible=================")
end

sossol=getvalue(sos)
hsol=getvalue(h)
fwsol=getvalue(fw)
aiwtsol=getvalue(aiwtUB)
cost=getobjectivevalue(m1st)
return hsol,sossol,fwsol,aiwtsol,cost
end
