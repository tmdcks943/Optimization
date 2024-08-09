#%%
import pandas as pd
import pandapower as pp
import numpy as np
from math import sqrt
import os
###############
#      README
# - pandapower에서는 name이 아닌 index 기반으로 참조함 ex) bus name 3이 net.bus의 두번째 행에 있을 경우(index 1) -> 1번 버스임
# - ieee 123 test feeder는 각 버스별 이름이 존재하나 index와 별개임
# - 본 코드에서는 각 버스별 이름을 인덱스로 설정하지 않았음
# - bus_to_index는 실제 ieee 123 feeder의 실제 bus 번호를 net.bus dataframe의 index로 대응하여 pandapower 에서 참조할 수 있는 이름으로 변경 
# - 실험용으로 모든 line에 switch를 추가하였음.
###############


def ieee123(with_der = True):
    previous_dir = os.getcwd()
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    line_data = pd.read_excel('IEEE123Node/line data.xls')
    line_config = pd.read_excel('IEEE123Node/line config.xls')

    net = pp.create_empty_network(name = 'IEEE123')
    for i in range(len(line_config)):
        pp.create_std_type(net, dict(line_config[line_config.columns[1:]].iloc[i]), name = line_config.name[i], element = 'line')

    bus_set = set(line_data[line_data.columns[0]]).union(set(line_data[line_data.columns[1]]))
    for i in bus_set:
        pp.create_bus(net, vn_kv = 4.16, name = f'{i}', id = i)
    ## external grid
    ext1 = pp.create_bus(net, vn_kv = 4.16, name = '150', id = 150)
    ext2 = pp.create_bus(net, vn_kv = 4.16, name = '233', id = 233)
#    ext3 = pp.create_bus(net, vn_kv = 4.16, name = '451', id = 451)
#    ext4 = pp.create_bus(net, vn_kv = 4.16, name = '195', id = 195)

    pp.create_ext_grid(net, ext1, vm_pu = 1.03, s_sc_max_mva = 3)
    pp.create_ext_grid(net, ext2, vm_pu = 1.03, s_sc_max_mva = 3)
#    pp.create_ext_grid(net, ext3, s_sc_max_mva = 5)
#    pp.create_ext_grid(net, ext4, s_sc_max_mva = 5)

    ## singlebus
    sb1 = pp.create_bus(net, vn_kv = 4.16, name = '610', id = 610)
    sb2 = pp.create_bus(net, vn_kv = 4.16, name = '350', id = 350)

    bus_to_index = pd.Series(net.bus.index, index = net.bus.name)         ## 이는 그림 상 나온 bus 이름 (int)을 받으면 그를 pandapower 상 index로 변환
    bus_to_index.index.name = None
    bus_to_index.index = bus_to_index.index.map(int)                

    ## capacitor
    
    capacitor_data = pd.read_excel('IEEE123Node/capacitor data.xls')
    for i in range(len(capacitor_data)):
        pp.create_shunt(net, bus_to_index.loc[capacitor_data[capacitor_data.columns[0]][i]], capacitor_data[capacitor_data.columns[1]][i]/1000.)

    ## line
    for i in range(len(line_data)): 
        pp.create_line(net, bus_to_index.loc[line_data[line_data.columns[0]][i]], bus_to_index.loc[line_data[line_data.columns[1]][i]], length_km = line_data[line_data.columns[2]][i]*0.0003048, std_type = line_data[line_data.columns[3]][i], name = f'{i+1}')
        pp.create_switch(net, bus_to_index.loc[line_data[line_data.columns[0]][i]], i, et = "l", closed = True) # 모든 line에 switch 연결
    pp.create_line(net, bus_to_index.loc[350], bus_to_index.loc[300], length_km = 0.0003, std_type = line_data[line_data.columns[3]][0], name = f'{i+2}')
    pp.create_line(net, bus_to_index.loc[61], bus_to_index.loc[610], length_km = 0.0003, std_type = line_data[line_data.columns[3]][0], name = f'{i+3}')



    load_data = pd.read_excel('IEEE123Node/spot loads data.xls')
    load_data['P'] = (load_data['P1'] + load_data['P2'] + load_data['P3'])/1000.
    load_data['Q'] = (load_data['Q1'] + load_data['Q2'] + load_data['Q3'])/1000.
    for i in range(len(load_data)):
        pp.create_load(net, bus = bus_to_index[load_data[load_data.columns[0]][i]], p_mw = load_data['P'][i], q_mvar = load_data['Q'][i], sn_mva = sqrt(load_data['P'][i]**2 + load_data['Q'][i]**2), name = f'{i+1}')
 

    switch_data = pd.read_excel('IEEE123Node/switch data.xls')
    for i in range(len(switch_data)):
        from_bus = switch_data[switch_data.columns[0]][i]
        to_bus = switch_data[switch_data.columns[1]][i]
        closed = switch_data[switch_data.columns[2]][i]
        if closed == "closed":
            closed = True
        else:
            closed = False
        pp.create_switch(net, bus_to_index.loc[from_bus], bus_to_index.loc[to_bus], et = "b", closed = closed)
    pp.create_switch(net, bus_to_index.loc[33], bus_to_index.loc[233], et = "b", closed = False) ## excel에 없음
    ##pp.create_switch(net, bus_to_index.loc[13], 12, et = "l", closed = True) ## excel에 없음
    pp.create_switch(net, bus_to_index.loc[39], bus_to_index.loc[66], et = "b", closed = False) ## excel에 없음
    pp.create_switch(net, bus_to_index.loc[17], bus_to_index.loc[96], et = "b", closed = False) ## excel에 없음
    ##pp.create_switch(net, bus_to_index.loc[76], 76, et = "l", closed = True) ## excel에 없음

    ## 참고 : 부하의 유효전력 합은 3.49 mw, 무효전력 합은 1.92 mvar

    sgen_data = pd.read_excel('IEEE123Node/sgen data.xls')

    if with_der is True:
        for i in range(len(sgen_data)):
            pp.create_sgen(net, bus_to_index[sgen_data[sgen_data.columns[0]][i]], p_mw=sgen_data[sgen_data.columns[1]][i], q_mvar=sgen_data[sgen_data.columns[2]][i], sn_mva=sgen_data[sgen_data.columns[3]][i], scaling=sgen_data[sgen_data.columns[4]][i], name=sgen_data[sgen_data.columns[5]][i], type=sgen_data[sgen_data.columns[6]][i])

    print("IEEE network is loaded")
    net.sgen.type[net.sgen.type == 'WT'] = 'WP'
    os.chdir(previous_dir)

    net.line.max_i_ka = net.line.max_i_ka*0.7
    return net
