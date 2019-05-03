Using the good cuts method:

Trial 1: all integer in first stage

Decompose the problem into two stage

first
    obj: cost
    var: x fw aiwt

second
    feasibility problem if aiwt is fixed, obj_sec will be fixed thru a series of continuous feasibility problem.

cutting idea: using nogoodcuts to cut aiwt out

Issue: slow

solution: search multiple region at the same for example ai=1/ai=0 separate a space into two.

result: 90% p=2 optimal solution is the UB with sum(fw)=20

        60% p=8 simulation in progress (have difficulty finish)

        20% p=16 simulation in progress (have difficulty finish)





Trial 2: only x and fw in first (x need to be binarilized or no-good-cuts can't be used here)
