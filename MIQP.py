import gurobipy as gp
from gurobipy import GRB
import pandapower as pp
import pandapower.networks as pn

import itertools
from Networks import ieee123

net = ieee123.ieee123() # P, Q를 일정값으로 하니까 풀림/ P, Q를 0으로 하면 loop 생김
#net = pn.case33bw()


#net = pn.create_cigre_network_mv(with_der=False)
Buses = list(net.bus.index)
SubstationBuses = list(net.ext_grid.bus)
NonsubstationBuses = list(net.bus.index)

for bus in SubstationBuses:
    NonsubstationBuses.remove(bus)

Lines = []
R = dict()
for idx in net.line.index:
    col = net.line.loc[idx]
    Lines.append((col.from_bus, col.to_bus))
    Lines.append((col.to_bus, col.from_bus))
    R[(col.from_bus, col.to_bus)] = col.r_ohm_per_km*col.length_km
    R[(col.to_bus, col.from_bus)] = col.r_ohm_per_km*col.length_km
# switch 만 연결된 line 도 line임

for idx in net.switch.index:
    col = net.switch.loc[idx]
    if col.et == 'b':
        Lines.append((col.bus, col.element))
        Lines.append((col.element, col.bus))
        R[(col.bus, col.element)] = 0.00001
        R[(col.element, col.bus)] = 0.00001

    if col.et == 't':
        from_bus = net.trafo.loc[net.switch.loc[idx].element].hv_bus
        to_bus = net.trafo.loc[net.switch.loc[idx].element].lv_bus
        
        if (from_bus, to_bus) not in Lines:
            Lines.append((from_bus, to_bus))
        if (to_bus, from_bus) not in Lines:
            Lines.append((to_bus, from_bus))
        
        R[(to_bus, from_bus)] = 0.00001
        R[(from_bus, to_bus)] = 0.00001
        
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
    P[bus] = 1
    Q[bus] = 1
    
for idx in net.load.index:
    bus = net.load.bus.loc[idx]
    p = net.load.p_mw.loc[idx]*net.load.scaling.loc[idx]
    q = net.load.q_mvar.loc[idx]*net.load.scaling.loc[idx]
    
    P[bus] += p
    Q[bus] += q
    
for idx in net.sgen.index:
    bus = net.sgen.bus.loc[idx]
    p = net.sgen.p_mw.loc[idx]*net.sgen.scaling.loc[idx]
    q = net.sgen.q_mvar.loc[idx]*net.sgen.scaling.loc[idx]
    
    P[bus] -= p
    Q[bus] -= q
    
TransferBuses = []
for bus in Buses:
    if P[bus] == 0 and Q[bus] == 0:
        TransferBuses.append(bus)

DistributedGeneratorBuses = []
for bus in Buses:
    if P[bus] < 0 or Q[bus] < 0:
        DistributedGeneratorBuses.append(bus)
        

model = gp.Model("MIQP_DNR") 

DistributedGeneratorBuses = []
for bus in Buses:
    if P[bus] < 0 or Q[bus] < 0:
        DistributedGeneratorBuses.append(bus)

p,q,x,y,k = dict(), dict(), dict(), dict(), dict()
for item in Lines:
    p[item] = model.addVar(name = 'p'+str(item))
    q[item] = model.addVar(name = 'q'+str(item))
    x[item] = model.addVar(vtype=GRB.BINARY, name = 'x'+str(item))
    #k[item] = model.addVar(name = 'k'+str(item))

#for bus in TransferBuses:
#    y[bus] = model.addVar(vtype=GRB.BINARY, name = 'y'+str(bus)) # eq 16
'''
K = dict()
for bus in Buses:
    if bus in DistributedGeneratorBuses:
        K[bus] = 1
    else:
        K[bus] = 0
'''

for bus in Buses: # eq5
    if not bus in SubstationBuses:
        model.addConstr(
            sum(
                p[(adjbus, bus)] - p[(bus, adjbus)] for adjbus in AdjacencyList[bus]
                ) == P[bus])

for bus in Buses: #eq6
    if not bus in SubstationBuses:
        model.addConstr(
            sum(
                q[(adjbus, bus)] - q[(bus, adjbus)] for adjbus in AdjacencyList[bus]
                ) == Q[bus])


#for bus in NonsubstationBuses:
#    model.addConstr(
#        -P[bus] - sum(x[(bus, adjbus)]*p[(bus, adjbus)] for adjbus in AdjacencyList[bus]) == 0
#    )
    
#for bus in NonsubstationBuses:
#    model.addConstr(
#        -Q[bus] - sum(x[(bus, adjbus)]*p[(bus, adjbus)] for adjbus in AdjacencyList[bus]) == 0
#    )
    
    
M = 10000000000
for line in Lines: 
    model.addConstr(p[line] >= 0) # eq 9 
    model.addConstr(p[line] <= M*x[line]) # eq 9
    model.addConstr(q[line] >= 0) # eq 10
    model.addConstr(q[line] <= M*x[line]) # eq 10

    #model.addConstr(x[line] >= 0) # eq 11


model.addConstr(
    sum(x[line] for line in Lines) == len(Buses) - len(SubstationBuses) #- sum(1 - y[bus] for bus in TransferBuses)
)
'''
for bus in TransferBuses:
    for adjbus in AdjacencyList[bus]:
        model.addConstr(
            x[(adjbus, bus)] <= y[bus]
        )
        model.addConstr(
            x[(bus, adjbus)] <= y[bus]
        )

    model.addConstr(
        sum(x[(bus, adjbus)] for adjbus in AdjacencyList[bus]) + sum(x[(adjbus, bus)] for adjbus in AdjacencyList[bus]) 
        >= 2*y[bus]
    )
'''
'''
for bus in Buses:
    model.addConstr(
        sum(k[adjbus, bus] for adjbus in AdjacencyList[bus])
        - sum(k[bus, adjbus] for adjbus in AdjacencyList[bus])
        == K[bus]
    )

for line in Lines:
    model.addConstr(
        k[line] <= len(DistributedGeneratorBuses)*x[line]
    )
    model.addConstr(
        -k[line] >= len(DistributedGeneratorBuses)*x[line]
    )
'''
model.setObjective(gp.quicksum(R[line]*(p[line]**2+q[line]**2) for line in Lines), GRB.MINIMIZE) # eq 4

model.update()

model.optimize()