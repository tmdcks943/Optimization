## MISOCP implementation of the following paper:
# Value of Distribution Network Reconfiguration in Presence of Renewable Energy Resources - II.B.

# Infeasiblity Scenarios
# 1. Large or small Sbase value
# 2. Large or small M values
# try relaxing constraint for 

#%%

import gurobipy as gp
from gurobipy import GRB
import pandapower as pp
import pandapower.networks as pn

import itertools
from Networks import ieee123

Vbase = 4.16e3 # 4.16kV
Sbase = 1e3 # 10MVA
Ibase = Sbase/Vbase #2.404e3 0.2404
Zbase = Vbase/Ibase #1.730e  1.730e4

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
    
    R[(col.from_bus, col.to_bus)] = col.r_ohm_per_km*col.length_km/Zbase
    R[(col.to_bus, col.from_bus)] =  col.r_ohm_per_km*col.length_km/Zbase
    
    X[(col.from_bus, col.to_bus)] = col.x_ohm_per_km*col.length_km/Zbase
    X[(col.to_bus, col.from_bus)] = col.x_ohm_per_km*col.length_km/Zbase
    

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
    P[bus] = 0.0001#*1000/Sbase
    Q[bus] = 0.0001#*1000/Sbase # 0.001 results in isolation

for idx in net.load.index:
    bus = net.load.bus.loc[idx]
    p = net.load.p_mw.loc[idx]*net.load.scaling.loc[idx]
    q = net.load.q_mvar.loc[idx]*net.load.scaling.loc[idx]
    
    P[bus] += p*1000/Sbase # kW /SBase
    Q[bus] += q*1000/Sbase
    
for idx in net.sgen.index:
    bus = net.sgen.bus.loc[idx]
    p = net.sgen.p_mw.loc[idx]*net.sgen.scaling.loc[idx]
    q = net.sgen.q_mvar.loc[idx]*net.sgen.scaling.loc[idx]
    
    P[bus] -= p*1000/Sbase
    Q[bus] -= q*1000/Sbase

model = gp.Model("MISOCP_DNR") 

#################
# Variable definition
p,q, l, z = dict(), dict(), dict(), dict()
vs = dict()
for line in Lines:
    p[line] = model.addVar(name = 'p'+str(line))
    q[line] = model.addVar(name = 'q'+str(line))

    l[line] = model.addVar(name = 'l'+str(line))


for line in Lines:
    z[line] = model.addVar(vtype=GRB.BINARY, name = 'z'+str(line))

for bus in NonsubstationBuses:
    vs[bus] = model.addVar(name = 'vs'+str(bus)) # v square
    #tvs[bus] = model.addVar(name = 'tvs'+str(bus)) # tilde v (defined for bus)

for bus in SubstationBuses: # eq 29
    vs[bus] = (4.16*1.03e3) ** 2/(Vbase ** 2) # 4.16kV 기준 1.03pu
    #tvs[bus]
    #tvs[bus] = model.addVar(name = 'tvs' + str(bus))

    P[bus] = model.addVar(name = 'P_F' + str(bus))
    Q[bus] = model.addVar(name = 'Q_F' + str(bus))
###################

# Objective function 
model.setObjective(sum(R[line]*l[line] for line in Lines), GRB.MINIMIZE) # eq9
#model.setObjective(sum(tp[line] + tq[line] for line in Lines), GRB.MINIMIZE) # eq9


################
# Radiality and substation flow 
for bus in SubstationBuses: 
    model.addConstr(
        sum(p[(bus, adjbus)] for adjbus in AdjacencyList[bus]) == P[bus]
    )
for bus in SubstationBuses: 
    model.addConstr(
        sum(q[(bus, adjbus)] for adjbus in AdjacencyList[bus]) == Q[bus]
    )
M = 100000 # i dont know why, but having larger M value results in infeasible situation.
# WhY? beacuse larger M value results in  M having large value
 
for line in Lines:  # eq 19
    model.addConstr(
        l[line] <= M*z[line] 
    )

model.addConstr( # Radiality eq 14 (modified)
    sum(z[line] for line in Lines) == len(Buses) - len(SubstationBuses)
)

for line in Switches:
    model.addConstr(
        z[line] + z[tuple(reversed(line))] <= 1
    )
for line in NoSwitches:
    model.addConstr(
        z[line] + z[tuple(reversed(line))] == 1
    )
##################
# SOCP flow equations
for bus in NonsubstationBuses: # eq 10
    model.addConstr(
        sum(R[(adjbus, bus)]*l[(adjbus, bus)] for adjbus in AdjacencyList[bus]) 
        ==
        -P[bus] +  sum(p[adjbus, bus] - p[bus, adjbus] for adjbus in AdjacencyList[bus])
    )

for bus in NonsubstationBuses: # eq 11
    model.addConstr(
        sum(X[(adjbus, bus)]*l[(adjbus, bus)] for adjbus in AdjacencyList[bus]) 
        ==
        -Q[bus] +  sum(q[adjbus, bus] - q[bus, adjbus] for adjbus in AdjacencyList[bus])
    )

M = 1000

for line in Lines: # eq 18 
    model.addConstr(
        (p[line]**2 + q[line]**2)
        <= vs[line[0]] * l[line] #line[1]
    )

for bus in NonsubstationBuses: # eq 27 
    for adjbus in AdjacencyList[bus]:
        model.addConstr(
            vs[bus]
            <= vs[adjbus] - 2*(R[(adjbus, bus)]*p[(adjbus, bus)] + X[(adjbus, bus)]*q[(adjbus, bus)]) + (R[(adjbus, bus)]**2 + X[(adjbus, bus)]**2)*l[(adjbus, bus)] + M*(1-z[adjbus, bus])
        )
for bus in NonsubstationBuses: # eq 28
    for adjbus in AdjacencyList[bus]:
        model.addConstr(
            vs[bus]
            >= vs[adjbus] - 2*(R[(adjbus, bus)]*p[(adjbus, bus)] + X[(adjbus, bus)]*q[(adjbus, bus)]) + (R[(adjbus, bus)]**2 + X[(adjbus, bus)]**2)*l[(adjbus, bus)]- M*(1-z[adjbus, bus])
        )

#model.setParam('MIPGap', 0.0001) # 1% gap
#model.setParam('FeasibilityTol', 1e-9)
#model.setParam('IntFeasTol', 1e-9)
model.setParam('Presolve', 0)

model.update()
model.optimize()


# %%
for i in vs:
    try:
        print(i, vs[i].X)
    except:
        pass
 #%%
# Flow from substation
vs[123] - 2*(R[(123,115)]*p[(123,115)].X + X[(123,115)]*q[(123,115)].X) - M*(1-z[123, 115].X)
#%%

#%%
# Flow from a bus connected to substation
vs[115].X - 2*(R[(115,0)]*p[(115,0)].X + X[(115,0)]*q[(115,0)].X)
#%%
M*(1-z[123, 115].X)
# %%
