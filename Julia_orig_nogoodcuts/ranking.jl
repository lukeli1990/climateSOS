function rankw(scenariow) #find the decreasing order of scenariow at each time and bus
# scenariow in B*T*sizew
B=14
T=3
sizew=20

orderscenw=zeros(Int,B,T,sizew)


        for i=1:B
		    for t=1:T
			orderscenw[i,t,:]=sortperm(scenariow[i,t,:],rev=true)
			end
		end
		return orderscenw
end
