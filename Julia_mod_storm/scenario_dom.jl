function scendom(scenariow)
    # Manually Changes: B sizew 

B=14
T=5
sizew=2000  #number of scenarios

domscen=zeros(Int,sizew,sizew) # for i!=j domscen(i,j)=1 ---> scen i >= scen j

for i=1:sizew-1
    for j=(i+1):sizew
        if minimum(scenariow[:,:,i]-scenariow[:,:,j])>=0
            domscen[i,j]=1
        end

        if minimum(scenariow[:,:,j]-scenariow[:,:,i])>=0
            domscen[j,i]=1
        end

    end
end

return domscen

end
