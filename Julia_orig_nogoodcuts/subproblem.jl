function subprob(a,scenind)
# a in B*T
# used in second stage to check fw
T=3

        for Tind=1:T
            println()
            println(" ===scenind: ",scenind," ===T: ",Tind)
            fwtemp=f2nd(a[:,Tind]) #simpler version of 2nd stage

				if fwtemp==0
				return 0
				end

		end
		return 1
end
