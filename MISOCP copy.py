## MISOCP implementation of the following paper:
# Convex models of distribution system reconfiguration (2012) - Section IV
################### eq 19, 7- 10, 
# Objective Function 
# eq 19 : the sum of the substation flows -> loss reduction objective : \sum_{i_\in B^F} p^F_i 
###################
# Constraints
###
# Radiality) eq 11 - 16 are replaced with N-S = L equation so that it can become feasible with supplying buses. z_{ij} = z_{ji} with 1 if connected and otherwise 0.
###
# Power flow constraints
# eq 21 - 29 : SOCP power flow equations
# \tilde{p}, \tilde{q} : total loss arose from supplying bus i ex) (j -> (loss)* -> i) sum of the losses
# eq 21, 22 : power balance constraint
# eq 23, 24 : calculate auxiliary voltage i from real voltage j -> utilized for calculating total loss \tilde{p}, \tilde{q} 
# eq 25, 26 : calculate total loss \tilde{p}, \tilde{q} 
# eq 27, 28 : calculate real voltage of bus i
# eq 29 : substation voltage
###
####################
##################
# Notice
# eq 25, 26 -> it seems wrong (it should be \sum_j {r_{ij} (p^2 + q^2)} <= \tilde{v}^2_i \tilde_p_i)
# replaced \tilde{p}_i and \tilde{q}_i with \tilde{p}_{ij} \tilde{q}_{ij} in eq 25, 26, and \tilde{p}_i = \sum_j \tilde{p}_{ij}
##################
#%%

Buses = [1,2,3,4]
NonsubstationBuses = [2,3,4]
SubstationBuses = [1,4]
Lines = [(1,2),(1,3),(2,1),(3,1), (2,3),(3,2), (2,4), (4,2)] 
Switches = [(1,2),(1,3),(2,1),(3,1), (2,3),(3,2), (2,4), (4,2)] 
NoSwitches = []

AdjacencyList = {1: [2,3], 2: [1,3,4], 3:[1,2], 4: [2]}

    # Parameters
P = {1 : 2.2, 2 : 1, 3 : 1, 4: 1} # P[1]은 사용되지 않음
Q = {1: 3, 2 : 2, 3 : 4, 4: 1} # Q[1]은 사용되지 않음
G = {(1,2) : 0, (2,1) : 0, (2,3) : 0, (3,2) : 0} # notused
B = {(1,2) : 10, (2,1) : 10, (2,3) : 10, (3,2) : 10} # notused
R = {(1,2) : 0.001, (2,1) : 0.001, (2,3) : 0.001, (3,2) : 0.001, (1,3) : 0.001, (3,1) :0.001, (2,4) : 0.001, (4,2) : 0.001}
X = {(1,2) : 0.001, (2,1) : 0.001, (2,3) : 0.001, (3,2) : 0.001, (1,3) : 0.001, (3,1) :0.001, (2,4) : 0.001, (4,2) : 0.001}



#%%

import gurobipy as gp
from gurobipy import GRB
import pandapower as pp
import pandapower.networks as pn

import itertools
from Networks import ieee123

Vbase = 4.16e3 # 4.16kV
Sbase = 1e7 # 10MVA
Ibase = Sbase/Vbase #2.404e3 2.404kA
Zbase = Vbase/Ibase #1.730

#net = pn.case33bw()
net = ieee123.ieee123()
#net = pn.case9()


Buses = list(net.bus.index)
SubstationBuses = list(net.ext_grid.bus)
NonsubstationBuses = list(net.bus.index)

for bus in SubstationBuses:
    NonsubstationBuses.remove(bus)

G, B = dict(), dict()
for bus in Buses:
    G[bus] = 0
    B[bus] = 0

Lines = []
R, X = dict(), dict()
for idx in net.line.index:
    col = net.line.loc[idx]
    Lines.append((col.from_bus, col.to_bus))
    Lines.append((col.to_bus, col.from_bus))
    
    R[(col.from_bus, col.to_bus)] = col.r_ohm_per_km*col.length_km #/Zbase
    R[(col.to_bus, col.from_bus)] =  col.r_ohm_per_km*col.length_km #/Zbase
    
    X[(col.from_bus, col.to_bus)] = col.x_ohm_per_km*col.length_km #/Zbase
    X[(col.to_bus, col.from_bus)] = col.x_ohm_per_km*col.length_km #/Zbase
    
# switch 만 연결된 line 도 line임

for idx in net.switch.index:
    col = net.switch.loc[idx]
    if col.et == 'b':
        Lines.append((col.bus, col.element))
        Lines.append((col.element, col.bus))
        R[(col.bus, col.element)] = 0.00001
        R[(col.element, col.bus)] = 0.00001
        
        X[(col.bus, col.element)] = 0.00001
        X[(col.element, col.bus)] = 0.00001

    if col.et == 't':
        from_bus = net.trafo.loc[net.switch.loc[idx].element].hv_bus
        to_bus = net.trafo.loc[net.switch.loc[idx].element].lv_bus
        
        if (from_bus, to_bus) not in Lines:
            Lines.append((from_bus, to_bus))
        if (to_bus, from_bus) not in Lines:
            Lines.append((to_bus, from_bus))
        
        R[(to_bus, from_bus)] = 0.00001
        R[(from_bus, to_bus)] = 0.00001

        X[(col.bus, col.element)] = 0.00001
        X[(col.element, col.bus)] = 0.00001
        
sw = dict()
for idx in net.switch.index:
    
    from_bus = net.switch.loc[idx].bus
    to_bus = None
    if net.switch.loc[idx].et == 'l':
        if net.line.loc[net.switch.loc[idx].element].to_bus != from_bus:
            to_bus = net.line.loc[net.switch.loc[idx].element].to_bus
        else:
            to_bus = net.line.loc[net.switch.loc[idx].element].from_bus
    if net.switch.loc[idx].et == 'b':
        to_bus = net.switch.loc[idx].element
    
    if net.switch.loc[idx].et == 't':
        to_bus = net.trafo.loc[net.switch.loc[idx].element].lv_bus
        
    sw[idx] = (from_bus, to_bus)
    
Switches = []
for switch in sw.values():
    if switch not in Switches:
        Switches.append(switch)
    if tuple(reversed(switch)) not in Switches:
        Switches.append(tuple(reversed(switch)))

from copy import deepcopy

Switches = []
for switch in sw.values():
    if switch not in Switches:
        Switches.append(switch)
    if tuple(reversed(switch)) not in Switches:
        Switches.append(tuple(reversed(switch)))

NoSwitches = deepcopy(Lines)

for switch in Switches:
    NoSwitches.remove(switch)


def build_adjacency_list(edges):
    adjacency_list = {}
    for node1, node2 in edges:
        if node1 not in adjacency_list:
            adjacency_list[node1] = []
        if node2 not in adjacency_list:
            adjacency_list[node2] = [] 
        #adjacency_list[node1].append(node2)
        adjacency_list[node2].append(node1)  # If the graph is undirected

    return adjacency_list
AdjacencyList = build_adjacency_list(Lines)

P = dict()
Q = dict()
for bus in Buses:
    P[bus] = 0.001#*1000/Sbase
    Q[bus] = 0.001#*1000/Sbase

for idx in net.load.index:
    bus = net.load.bus.loc[idx]
    p = net.load.p_mw.loc[idx]*net.load.scaling.loc[idx]
    q = net.load.q_mvar.loc[idx]*net.load.scaling.loc[idx]
    
    P[bus] += p#*1000/Sbase
    Q[bus] += q#*1000/Sbase
    
for idx in net.sgen.index:
    bus = net.sgen.bus.loc[idx]
    p = net.sgen.p_mw.loc[idx]*net.sgen.scaling.loc[idx]
    q = net.sgen.q_mvar.loc[idx]*net.sgen.scaling.loc[idx]
    
    P[bus] -= p*5 #*1000/Sbase
    Q[bus] -= q*5 #*1000/Sbase

model = gp.Model("MISOCP_DNR") 

#################
# Variable definition
p,q, tp, tq, z = dict(), dict(), dict(), dict(), dict()
vs, tvs = dict(), dict()
for line in Lines:
    p[line] = model.addVar(name = 'p'+str(line))
    q[line] = model.addVar(name = 'q'+str(line))

    tp[line] = model.addVar(name = 'tp'+str(line)) # tilde p (defined for line)
    tq[line] = model.addVar(name = 'tq'+str(line)) # tilde q (defined for line)


for line in Lines:
    z[line] = model.addVar(vtype=GRB.BINARY, name = 'z'+str(line))

for bus in NonsubstationBuses:
    vs[bus] = model.addVar(name = 'vs'+str(bus)) # v square
    tvs[bus] = model.addVar(name = 'tvs'+str(bus)) # tilde v (defined for bus)

for bus in SubstationBuses: # eq 29
    vs[bus] = 1.03 ** 2 # 4.16 ** 2
    #tvs[bus] = 1.03 ** 2
    tvs[bus] = model.addVar(name = 'tvs' + str(bus))

    P[bus] = model.addVar(name = 'P_F' + str(bus))
    Q[bus] = model.addVar(name = 'Q_F' + str(bus))
###################

# Objective function : minimize substation flow
model.setObjective(sum(tp[line] + tq[line] for line in Lines), GRB.MINIMIZE) # eq9
#model.setObjective(sum(tp[line] + tq[line] for line in Lines), GRB.MINIMIZE) # eq9


################
# Radiality and substation flow (eq 7- 16)
for bus in SubstationBuses: # eq 7
    model.addConstr(
        sum(p[(bus, adjbus)] for adjbus in AdjacencyList[bus]) == P[bus]
    )
for bus in SubstationBuses: # eq 8
    model.addConstr(
        sum(q[(bus, adjbus)] for adjbus in AdjacencyList[bus]) == Q[bus]
    )
M = 1000
 
for line in Lines: # eq 9
    model.addConstr(
        p[line] >= 0 
    )
    model.addConstr(
        p[line] <= M*z[line]
    )
for line in Switches: # eq 10
    model.addConstr(
        q[line] >= 0 
    )
    model.addConstr(
        q[line] <= M*z[line]
    )

for line in Lines: # eq 9
    model.addConstr(
        tp[line] >= 0 
    )
    model.addConstr(
        tp[line] <= M*z[line]
    )
for line in Lines: # eq 10
    model.addConstr(
        tq[line] >= 0 
    )
    model.addConstr(
        tq[line] <= M*z[line]
    )


model.addConstr( # Radiality
    sum(z[line] for line in Lines) == len(Buses) - len(SubstationBuses)
)

for line in Switches:
    model.addConstr(
        z[line] + z[tuple(reversed(line))] <= 1
    )
##################
# SOCP flow equations
for bus in NonsubstationBuses: # eq 21
    model.addConstr(
        sum(tp[(adjbus, bus)] for adjbus in AdjacencyList[bus]) 
        ==
        -P[bus] +  sum(p[adjbus, bus] - p[bus, adjbus] for adjbus in AdjacencyList[bus])
    )

for bus in NonsubstationBuses: # eq 22
    model.addConstr(
        sum(tq[(adjbus, bus)] for adjbus in AdjacencyList[bus]) 
        ==
        -Q[bus] +  sum(q[adjbus, bus] - q[bus, adjbus] for adjbus in AdjacencyList[bus])
    )

M = 1000
for bus in Buses: # eq 23
    for adjbus in AdjacencyList[bus]:
        model.addConstr(
            tvs[adjbus] 
            <= vs[bus] + M*(1-z[bus, adjbus])
        )

for bus in Buses: # eq 24
    for adjbus in AdjacencyList[bus]:
        model.addConstr(
            tvs[adjbus] 
            >= vs[bus] - M*(1-z[bus, adjbus])
        )
for line in Lines: # eq 25 (modified)
    model.addConstr(
        R[line]*(p[line]**2 + q[line]**2)
        <= tvs[line[0]] * tp[line] #line[1]
    )

for line in Lines: # eq 26 (modified)
    model.addConstr(
        X[line]*(p[line]**2 + q[line]**2)
        <= tvs[line[0]] * tq[line]
    )
for bus in NonsubstationBuses: # eq 27 #Nonsub
    for adjbus in AdjacencyList[bus]:
        model.addConstr(
            vs[bus]
            <= vs[adjbus] - 2*(R[(adjbus, bus)]*p[(adjbus, bus)] + X[(adjbus, bus)]*q[(adjbus, bus)]) + M*(1-z[adjbus, bus])
        )
for bus in NonsubstationBuses: # eq 28
    for adjbus in AdjacencyList[bus]:
        model.addConstr(
            vs[bus]
            >= vs[adjbus] - 2*(R[(adjbus, bus)]*p[(adjbus, bus)] + X[(adjbus, bus)]*q[(adjbus, bus)]) - M*(1-z[adjbus, bus])
        )
#for bus in Buses: # eq 27 #Nonsub
#    for adjbus in AdjacencyList[bus]:
#        model.addConstr(
#            vs[adjbus]
#            <= vs[bus] - 2*(R[(adjbus, bus)]*p[(bus, adjbus)] + X[(adjbus, bus)]*q[(bus, adjbus)]) + M*(1-z[bus, adjbus])
#        )
#for bus in Buses: # eq 28
#    for adjbus in AdjacencyList[bus]:
#        model.addConstr(
#            vs[adjbus]
#            >= vs[bus] - 2*(R[(adjbus, bus)]*p[(bus, adjbus)] + X[(adjbus, bus)]*q[(bus, adjbus)]) - M*(1-z[bus, adjbus])
#        )
#model.setParam('MIPGap', 0.0001) # 1% gap
#model.setParam('FeasibilityTol', 1e-9)
#model.setParam('IntFeasTol', 1e-9)
#model.setParam('Presolve', 0)

model.update()
model.optimize()


# %%
for i in vs:
    print(i, vs[i].X)
 #%%
# Flow from substation
vs[123] - 2*(R[(123,115)]*p[(123,115)].X + X[(123,115)]*q[(123,115)].X)

#%%
# Flow from a bus connected to substation
vs[115].X - 2*(R[(115,0)]*p[(115,0)].X + X[(115,0)]*q[(115,0)].X)
#%%
M*(1-z[adjbus, bus].X)
# %%
#############################
# Regarding Infeasiblity
# Infeasible situation occured due to the following reasons.
# Reason 1) z being slightly nonzero and M being very big -> eq 27 and 28 did not bind voltages
# Solution) Set appropriate M value. You can relax other constraints to see M(1-z) values.

# Reason 2) Also check if P, Q, R, X is too large so that voltage deviates so badly.
# Large P (power injection) -> large p (branch power) 
# Large P and large R -> large Rp -> eq 27 voltage drop
# Solution) Set appropriate P, Q, R, X value. (Decrease P, Q, R, X)

# Experiment Result
# It is fast but it can't consider transfer node. 0.1 minimal bus injection results in one bus votlage being larger than 220, and 0.01 minimal bus injection results in ltos of bus voltage being larger than 220.
##############################

##############################
# Problems of the original paper
# 1. Wrong notation on (25, 26)
# The inequality is a function of j, which results in multiple constraint on single i
# Rather, it should be \tilde{p}_{ji}, and \tilde{p}_i in (21) should be replaced with sum_j \tilde{p}_ji

# 2. Not including q in objective function
# Not including q in minimzation results in large q tilde value due to lack of minimizing effect on q
# larger q -> small v (27, 28) -> small v tilde (23, 24) -> small p tilde (25)
#############################