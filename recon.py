import gurobipy as gp
from gurobipy import GRB
import itertools
import numpy as np
import tqdm as tqdm
from copy import deepcopy


class MIQP():
    def __init__(self, env, name = 'MIQP'):
        self.env = env
        self.net = env.net
        self.name = 'MIQP'

    def act(self):
        return
    def SolveProblem(self):
        model = gp.Model("MIQP_DNR") 

        p,q,y,z = dict(), dict(), dict(), dict()
        for item in self.Lines:
            p[item] = model.addVar(name = 'p'+str(item))
            q[item] = model.addVar(name = 'q'+str(item))
            z[item] = model.addVar(name = 'z'+str(item))

        for line in self.Switches:
            y[line] = model.addVar(vtype=GRB.BINARY,name = 'y'+str(line)) # eq 16

            
        # Equality Constraint
        for bus in self.Buses: # eq5
            if not bus in self.SubstationBuses:
                model.addConstr(
                    sum(
                        p[(adjbus, bus)] - p[(bus, adjbus)] for adjbus in self.AdjacencyList[bus]
                        ) == self.P[bus])

        for bus in self.Buses: #eq6
            if not bus in self.SubstationBuses:
                model.addConstr(
                    sum(
                        q[(adjbus, bus)] - q[(bus, adjbus)] for adjbus in self.AdjacencyList[bus]
                        ) == self.Q[bus])
        # Nonequality Constraint
        M = 1000000
        for line in self.Lines: 
            model.addConstr(p[line] >= 0) # eq 9 
            model.addConstr(p[line] <= M*z[line]) # eq 9
            model.addConstr(q[line] >= 0) # eq 10
            model.addConstr(q[line] <= M*z[line]) # eq 10

            model.addConstr(z[line] >= 0) # eq 11
            if line[1] in self.SubstationBuses:
                model.addConstr(z[line] == 0) # eq 12 # 이거의 의미 : 역송 방향 (bus -> substation)은 없다 이거잖아

        for line in self.NoSwitches:
            model.addConstr(z[line] + z[tuple(reversed(line))] == 1) # eq 13

        for line in self.Switches:
            model.addConstr(z[line] + z[tuple(reversed(line))] == y[line]) # eq 14
            
        for bus in self.Buses:
            if not bus in self.SubstationBuses:
                model.addConstr(
                    sum(
                        z[(adjbus, bus)] for adjbus in self.AdjacencyList[bus]
                    ) == 1 
                ) # eq 15

        model.setObjective(gp.quicksum(self.R[line]*(p[line]**2+q[line]**2) for line in self.Lines), GRB.MINIMIZE) # eq 4
        model.update()
        model.optimize()

        if model.status == GRB.OPTIMAL:
            action = np.zeros(shape = (len(self.net.switch)))
            for idx, switch in enumerate(self.Switches):
                action[idx//2] = y[switch].X
            action = np.round(action).astype(int)
        else:
            action = np.zeros(shape = (len(self.net.switch)))
            action = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1]
        return action

    def evaluation(self, eval_len = 288,  use_tqdm = True): # test trajectory에서 테스트       
        self.time_elapsed = 0
        for i in tqdm(range(eval_len)):
            self.env.iter = i
            action = self.act()
            _, reward = self.env.step(action, save_result = True, maintain_data = False)
        now = datetime.now()
        self.env.ResetSeed(now.microsecond)
        self.env.SaveResult.save_result(name = self.name + now.strftime('%Y%m%d%H%M%S'))
        return self.env.SaveResult.result()
    
    def VariableInitialization(self):
        net = self.net

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


        self.Buses = Buses
        self.SubstationBuses = SubstationBuses
        self.NonsubstationBuses = NonsubstationBuses
        self.Lines = Lines
        self.R = R
        self.Switches = Switches
        self.NoSwitches = NoSwitches
        self.AdjacencyList = AdjacencyList

    def LoadPQ(self):      
        Buses = self.Buses
        net = self.net

        P = dict()
        Q = dict()
        for bus in Buses:
            P[bus] = 0
            Q[bus] = 0
            
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
        self.P, self.Q = P, Q