function scenmanipulate(scenariow,conflvl,critbus)
    # Manually Changes: B sizew 

B=14
T=5
sizew=2000  #number of scenarios
thresind0=ceil(Int, conflvl*sizew)

rankscen=zeros(Int,B,T,sizew)
scenariosort=zeros(B,T,sizew)
thresind=zeros(Int,B,T)

for i=1:B
    for t=1:T
        scenariosort[i,t,:]=sort(scenariow[i,t,:]) #ascend
        rankscen[i,t,:]=sortperm(scenariow[i,t,:]) #ascend ind
            if critbus[i,t]==1 && scenariosort[i,t,thresind0]>0
                thresind[i,t]=thresind0
            else
               for s=1:sizew
                   if scenariosort[i,t,s]==0
                       thresind[i,t]=s
                   end
               end
            end
        end
    end

if minimum(thresind)==0
    println("there exists minimum scenario that is not zero")
end


return thresind,rankscen,scenariosort


end
