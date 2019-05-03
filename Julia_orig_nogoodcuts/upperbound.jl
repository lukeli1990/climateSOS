function upperboundx(scenariow,ei,di)#calculate the upperbound for each bus
# scenariow in T*sizew;ei di in scalar
# used in cutting to set
T=3
sizew=20
UB=0

        for Tind=1:T
		    for scenind=1:sizew
			UB=max(ceil((scenariow[Tind,scenind]-ei)/di),UB)
			end
		end
		return UB
end
