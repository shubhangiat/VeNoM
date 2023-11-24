from os import listdir

def adeg(e, v):
    return (e*2.0/v)

deginfo = open('deg_split_info', 'w')

for i in range(1, 21):
    with open('exp_l150_dense/q_v7/q_v7_edges_' + str(i)) as f:
        raw=f.readlines()
    v = 7
    e = len(raw)
    deg = adeg(e, v)
    deg = round(deg, 2)
    # if deg < 2:
    #     d = 'deg12'
    # elif deg < 3:
    #     d = 'deg23'
    # elif deg < 4:
    #     d = 'deg34'
    # elif deg < 5:
    #     d = 'deg45'
    # elif deg < 6:
    #     d = 'deg56'
    # elif deg < 7:
    #     d = 'deg67'
    deginfo.write('exp_l150_dense/q_v7/q_v7_edges_' + str(i) + '\t' + str(deg) + '\n')

for i in range(1, 21):
    with open('exp_l150_dense/q_v9/q_v9_edges_' + str(i)) as f:
        raw=f.readlines()
    v = 9
    e = len(raw)
    deg = adeg(e, v)
    deg = round(deg, 2)
    # if deg < 2:
    #     d = 'deg12'
    # elif deg < 3:
    #     d = 'deg23'
    # elif deg < 4:
    #     d = 'deg34'
    # elif deg < 5:
    #     d = 'deg45'
    # elif deg < 6:
    #     d = 'deg56'
    # elif deg < 7:
    #     d = 'deg67'
    # elif deg < 8:
    #     d = 'deg78'
    # elif deg < 9:
    #     d = 'deg89'
    deginfo.write('exp_l150_dense/q_v9/q_v9_edges_' + str(i) + '\t' + str(deg) + '\n')

    
for i in range(1, 21):
    with open('exp_l150_sparse/q_v7/q_v7_edges_' + str(i)) as f:
        raw=f.readlines()
    v = 7
    e = len(raw)
    deg = adeg(e, v)
    deg = round(deg, 2)
    # if deg < 2:
    #     d = 'deg12'
    # elif deg < 3:
    #     d = 'deg23'
    # elif deg < 4:
    #     d = 'deg34'
    # elif deg < 5:
    #     d = 'deg45'
    # elif deg < 6:
    #     d = 'deg56'
    # elif deg < 7:
    #     d = 'deg67'
    deginfo.write('exp_l150_sparse/q_v7/q_v7_edges_' + str(i) + '\t' + str(deg) + '\n')

for i in range(1, 21):
    with open('exp_l150_sparse/q_v9/q_v9_edges_' + str(i)) as f:
        raw=f.readlines()
    v = 9
    e = len(raw)
    deg = adeg(e, v)
    deg = round(deg, 2)
    # if deg < 2:
    #     d = 'deg12'
    # elif deg < 3:
    #     d = 'deg23'
    # elif deg < 4:
    #     d = 'deg34'
    # elif deg < 5:
    #     d = 'deg45'
    # elif deg < 6:
    #     d = 'deg56'
    # elif deg < 7:
    #     d = 'deg67'
    # elif deg < 8:
    #     d = 'deg78'
    # elif deg < 9:
    #     d = 'deg89'
    deginfo.write('exp_l150_sparse/q_v9/q_v9_edges_' + str(i) + '\t' + str(deg) + '\n')

    
deginfo.close()