println()
println()
println()
println()
using MatpowerCases
using JLD
include("find_tree.jl")
# presolve order 2
# redistribute load based on sealevel rize; Type II in the documentation
# Manually Changes: B sizew jldname
     T=5
     B=14
     sizew=100

     submind=zeros(B,T,sizew)
     loadrd=zeros(B,T,sizew)

     scendata=load("14scenario100+sea.jld")

     mpc = loadcase("case14")

    #identify the load index and number
     loadval0=mpc["bus"][:,3] #B
     loadnum=0
     for i=1:B
         if loadval0[i]>0
          loadnum=loadnum+1
         end
     end
     loadind=zeros(Int, 1,loadnum)
     tempind=1
     for i=1:B
         if loadval0[i]>0
          loadind[tempind]=i
          tempind=tempind+1
         end
     end
     if loadnum != tempind-1
          println("load count mistake")
     end


     ei=scendata["ei"]'*ones(1,T)#B*T
     SL=scendata["SLscen"] #B*T*sizew

    for s=1:sizew
         for t=1:T
              for i=1:B
                    if ei[i,t]<SL[i,t,s]
                         submind[i,t,s]=1
                    end
               end
               if minimum(submind[:,t,s])==1
                    println("all buses are submerged: ",t,";",s)
               end
          end
     end


     #generate tree and branch for particular bus
     treenet=zeros(Int,B,B,B)
     for i=1:B
          (treenet[i,:,:])=findtree(i)
     end



    #loadredistribution
    for s=1:sizew
         for t=1:T
              loadrd[:,t,s]=loadval0 #initialization on load values
              if maximum(submind[:,t,s])>0 #there are submerged bus
                   for i=1:B
                        if submind[i,t,s]==1 && loadrd[i,t,s]>0 #find the submerged load bus
                             busexist=(1-submind[:,t,s])*ones(1,B) # 1 exsit 0 subumerge
                             treexist=busexist.*treenet[i,:,:]
                             for level=2:B
                                  if maximum(treexist[:,level])>0
                                       loadrd[:,t,s]=loadrd[:,t,s]+treexist[:,level]*loadrd[i,t,s]/sum(treexist[:,level])
                                       loadrd[i,t,s]=0
                                       treexist=zeros(Int,B,B)
                                  end
                             end
                             if maximum(treexist)>0
                                  println()
                                  println()
                                  println("ERROR on load redistribution")
                             end

                        end
                   end
              end
         end
    end

   #Unfinish
   #filter out the sea level failed cases assumption A2 in Luedtke
