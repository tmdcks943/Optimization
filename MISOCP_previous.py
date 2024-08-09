# MICP implementation of the folloing paper:
# Value of Distribution Network Reconfiguration in Presence of Renewable Energy Resources (2016)


import gurobipy as gp
from gurobipy import GRB
import pandapower as pp
import pandapower.networks as pn

import itertools
from Networks import ieee123

#net = pn.case33bw()
net = ieee123.ieee123()
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
    
    R[(col.from_bus, col.to_bus)] = col.r_ohm_per_km*col.length_km
    R[(col.to_bus, col.from_bus)] = col.r_ohm_per_km*col.length_km
    
    X[(col.from_bus, col.to_bus)] = col.x_ohm_per_km*col.length_km
    X[(col.to_bus, col.from_bus)] = col.x_ohm_per_km*col.length_km
    
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
    P[bus] = -1
    Q[bus] = -1
    
for idx in net.load.index:
    bus = net.load.bus.loc[idx]
    p = net.load.p_mw.loc[idx]*net.load.scaling.loc[idx]
    q = net.load.q_mvar.loc[idx]*net.load.scaling.loc[idx]
    
    P[bus] -= p
    Q[bus] -= q
    
for idx in net.sgen.index:
    bus = net.sgen.bus.loc[idx]
    p = net.sgen.p_mw.loc[idx]*net.sgen.scaling.loc[idx]
    q = net.sgen.q_mvar.loc[idx]*net.sgen.scaling.loc[idx]
    
    P[bus] -= 0#p test용
    Q[bus] -= 0#q
    
    

model = gp.Model("MISOCP_DNR") 

p,q,l,a = dict(), dict(), dict(), dict()
v = dict()
for line in Lines:
    p[line] = model.addVar(name = 'p'+str(line))
    q[line] = model.addVar(name = 'q'+str(line))
    l[line] = model.addVar(name = 'l'+str(line))
    a[line] = model.addVar(vtype=GRB.BINARY, name = 'a'+str(line))

for bus in Buses:
    v[bus] = model.addVar(name = 'v'+str(bus))

for bus in NonsubstationBuses: # eq10
    model.addConstr(
        sum(
            p[(bus, adjbus)] for adjbus in AdjacencyList[bus]
            ) - 
        sum(
            p[(adjbus, bus)] - R[(adjbus, bus)]*l[(adjbus, bus)] for adjbus in AdjacencyList[bus]
            )
        + G[bus]*v[bus]
        == -P[bus])

for bus in NonsubstationBuses: # eq11
    model.addConstr(
        sum(
            q[(bus, adjbus)] for adjbus in AdjacencyList[bus]
            ) - 
        sum(
            q[(adjbus, bus)] - X[(adjbus, bus)]*l[(adjbus, bus)] for adjbus in AdjacencyList[bus]
            )
        + B[bus]*v[bus]
        == -Q[bus])

M = 100000
for line in Switches: # eq12
    bus1, bus2 = line[0], line[1]
    model.addConstr(
        M*(1-a[line]) + v[bus1] -  2*(R[line]*p[line] + X[line]*q[line]) + (R[line]**2 + X[line]**2)*l[line]
        >= v[bus2]
    )
    
for line in Switches: # eq13
    bus1, bus2 = line[0], line[1]
    model.addConstr(
        -M*(1-a[line]) + v[bus1] -  2*(R[line]*p[line] + X[line]*q[line]) + (R[line]**2 + X[line]**2)*l[line]
        <= v[bus2]
    )
     
model.addConstr( # eq 14
    sum(
        a[line] for line in Switches
    ) == len(Buses) - len(SubstationBuses)
)


for line in Switches: # eq 18
    model.addConstr(
        (p[line]**2 + q[line]**2)
        <= l[line]*v[line[0]]
    )

for line in Switches: # eq 19
    model.addConstr(
        a[line]*M
        >= l[line]
    )
    
for bus in SubstationBuses: # eq 20
    model.addConstr(
        v[bus] == 1.03
    )

for bus in SubstationBuses: # eq 21
    model.addConstr(
        v[bus] <= 2
    )
    
for bus in NonsubstationBuses: # eq 21
    model.addConstr(
        v[bus] >= 0.4
    )

model.setObjective(sum(R[line]*l[line] for line in Switches), GRB.MINIMIZE) # eq9
model.update()
model.optimize()