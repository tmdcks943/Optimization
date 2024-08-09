import gurobipy as gp
from gurobipy import GRB
import itertools

try:
    # Radial power flow CP (Radial Distribution Load Flow Using Conic Programming, 2006)
    # decision variables
    # u = V^2/sqrt(2) (bus에 대한 variable)
    # R = V_i V_j cos \theta, I = v_i V_j sin \theta (line 에 대한 variable) /  single substation 에서 bus 개수 - 1 개의 line을 가짐
    
    # Create a new model
    model = gp.Model("powerflow_CP") 

    # Define the sets of generators and buses (nodes)
    buses = ['B1', 'B2', 'B3']
    NonsubstationBuses = ['B2', 'B3']
    SubstationBuses = ['B1']
    lines = ['L1', 'L2'] 

    # Incident lines
    IncidentLines = {'B2' : ['L1', 'L2'], 'B3' : ['L2']}
    
    # AdjacencyList
    AdjacencyList = {'B1' : ['B2'], 'B2' : ['B1', 'B3'], 'B3' : ['B2']}  # Adjacency List
        
    # Parameters
    P = {'B2' : 1, 'B3' : 1}
    Q = {'B2' : 0, 'B3' : 0}
    G = {('B1', 'B2') : 0, ('B2', 'B1') : 0, ('B2', 'B3') : 0, ('B3', 'B2') : 0}
    B = {('B1', 'B2') : 10, ('B2', 'B1') : 10, ('B2', 'B3') : 10, ('B3', 'B2') : 10}
    

    # Decision Variable ( Substation bus는 u 가 그냥 const이고, NonsubstationBus는 Variable - *addVars가 리턴하는건 Gurobi Tupledict로, 그 일부를 const로 지정 불가능하여 일반 dict로 initialize)
    u = dict()
    for SubstationBus in SubstationBuses:
        u[SubstationBus] = 1/2**(1/2) # V1 = 1 (eq 10)
    for NonsubstationBus in NonsubstationBuses:
        u[NonsubstationBus] = model.addVar(name = 'u(' + NonsubstationBus +')')
    
    
    R, I = dict(), dict()
    for Bus in buses:
        for AdjacentBus in AdjacencyList[Bus]:
            if (AdjacentBus, Bus) in R.keys():
                R[(Bus, AdjacentBus)] = R[(AdjacentBus, Bus)] 
                I[(Bus, AdjacentBus)] = -I[(AdjacentBus, Bus)] 
            else:
                R[(Bus, AdjacentBus)] = model.addVar(name = 'R(' + Bus + ','+AdjacentBus+')')
                I[(Bus, AdjacentBus)] = model.addVar(name = 'I(' + Bus + ','+AdjacentBus+')')

    # Equality Constraints (Power flow equation)
    for Bus in NonsubstationBuses:
        model.addConstr(P[Bus] == -(2**(1/2) * u[Bus] * sum(G[(Bus, AdjacentBus)] for AdjacentBus in AdjacencyList[Bus])) + sum(G[(Bus, AdjacentBus)]*R[(Bus, AdjacentBus)] - B[(Bus, AdjacentBus)]*I[(Bus, AdjacentBus)] for AdjacentBus in AdjacencyList[Bus]), "eq6") # eq 6
        model.addConstr(Q[Bus] == -(2**(1/2) * u[Bus] * sum(B[(Bus, AdjacentBus)] for AdjacentBus in AdjacencyList[Bus])) + sum(B[(Bus, AdjacentBus)]*R[(Bus, AdjacentBus)] + G[(Bus, AdjacentBus)]*I[(Bus, AdjacentBus)] for AdjacentBus in AdjacencyList[Bus]), "eq6") # eq 6
    
    # Nonequality Constraints
    checker = dict()
    for Bus in buses:
        for AdjacentBus in AdjacencyList[Bus]:            
            if not (AdjacentBus,Bus) in checker.keys():
                checker[(Bus, AdjacentBus)] = True
                model.addConstr(2*u[Bus]*u[AdjacentBus] >= R[(Bus, AdjacentBus)]**2 + I[(Bus, AdjacentBus)]**2, 'eq9')         # eq 9    
                model.addConstr(R[(Bus, AdjacentBus)] >= 0, 'eq11')

    for Bus in NonsubstationBuses:
        model.addConstr(u[Bus] >=  0, 'eq10') # eq 10
    
    # Objective : Maximize R
    model.setObjective(gp.quicksum(R[(Bus, AdjacentBus)] for Bus in buses for AdjacentBus in AdjacencyList[Bus]), GRB.MAXIMIZE)
    
    # Optimize the model
    model.optimize()
    print(model.getVars())
    print(model.getConstrs())

    for v in model.getVars():
        print('%s %g' % (v.varName, v.x))
        
except gp.GurobiError as e:
    print(f'Error code {e.errno}: {e}')

except AttributeError:
    print('Encountered an attribute error')
