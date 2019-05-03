modified from Luedtke's method by moving all integer variables into first stage and combine as master problem. No decomposition
except check for the feasibilities.

Decompose the problem into two stage

first 
    obj: cost 
    var: x fw aiwt

second
    feasibility problem if aiwt is fixed, obj_sec will be fixed thru a series of continuous feasibility problem. 

cutting idea: modified from Luedtke's method; see documentation 

Issue: in progress