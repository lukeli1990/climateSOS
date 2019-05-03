function subprob(a,scenind)
# a in B*T
# used in second stage to check fw
T=3

        for Tind=1:T

			
            fwtemp=f2nd(a[:,Tind]) #simpler version of 2nd stage
				#println(" ===scenind: ",scenind," ===T: ",Tind)

				if fwtemp==0
                    println(" =====scenind: ",scenind," fail the current design at t=",Tind,"=========")
                    println()
				return 0,Tind
				end

		end
        println(" =====scenind: ",scenind," pass the current design =========")
        println()
		return 1,T
end
