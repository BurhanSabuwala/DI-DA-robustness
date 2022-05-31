import numpy as np
from scipy.stats import zscore


def get_core_network(name):
    """
    This function reads the topo files and fetches names of the nodes, core nodes (nodes other than signalling nodes and output nodes), signalling nodes (no incoming edges), output nodes (no outgoing edges) and indexing dictionary 
    """
    content = open(name+'.topo').read().split('\n')
    content = [content[i].split() for i in range(len(content))]
    content.remove(content[0])
    name_nodes =[]
    ind = 0
    index_nodes = {}
    for i in range(len(content)):
        if len(content[i])!=3:
            continue
        if content[i][0] not in name_nodes:
            name_nodes.append(content[i][0])
            index_nodes[content[i][0]] = ind
            ind+=1
        if content[i][1] not in name_nodes:
            name_nodes.append(content[i][1])
            index_nodes[content[i][1]] = ind
            ind = ind + 1
    nodes = len(name_nodes)
    sig_nodes = [i for i in name_nodes]
    out_nodes = [i for i in name_nodes]
    """ Identifying Output and signalling nodes"""
    for i in range(len(content)):
        if len(content[i])!=3:
            continue
        if content[i][0] in out_nodes:
            out_nodes.remove(content[i][0])
            
        if content[i][1] in sig_nodes:
            sig_nodes.remove(content[i][1])
    core_nodes = []
    for i in range(len(name_nodes)):
        if name_nodes[i] not in out_nodes and name_nodes[i] not in sig_nodes:
            core_nodes.append(name_nodes[i])
    name_nodes = core_nodes + sig_nodes + out_nodes
    return name_nodes, core_nodes, sig_nodes, out_nodes, index_nodes


def ss_zscore(name_of_topo, core_nodes, index_nodes):
    """
    This function is used to analyse the RACIPE output. It 
    normalizes the steady states values that RACIPE outputs 
    with respect to production rate and degradation rates of 
    the incoming edges. This is followed by z-score calculation 
    on the normalized values of the steady states.
    """
    print(name_of_topo)
    nodes= len(index_nodes)
    params = np.loadtxt(str(name_of_topo)+"_parameters.dat")
    ss_v = [params[i][1] for i in range(len(params))]
    ss_v = list(np.sort(list(set(ss_v))))
    ss = [np.loadtxt(str(name_of_topo)+"_solution_"+str(int(ssv))+".dat") for ssv in ss_v]
    nor =[[1 for i in range(len(params))]for i in range(nodes)]
    f = open(str(name_of_topo) + '.prs', 'r')
    prs = f.read(); f.close();
    prs = prs.split('\n')
    prs.remove(prs[0])
    prs = [prs[i].split('\t') for i in range(len(prs))]
    prs = [prs[i][0] for i in range(len(prs))]
    prs = [prs[i].split('To') for i in range(len(prs))]
    
    """Creating Normalization Matrix"""
    index_fav = [2];
    zs_fav = [[] for i in range(nodes)]
    for i in range(len(params)):
        j=2
        while(j<len(params[i])):
            if j<(nodes+2):
                if 'Prod' in prs[j-2][0]:
                    nor[j-2][i] = nor[j-2][i]*params[i][j]
                #print(j, prs[j-2])
                j=j+1
                
            elif j<(2*nodes+2):
                ni = j-nodes-2
                if 'Deg' in prs[j-2][0]:
                    nor[ni][i] = nor[ni][i]/params[i][j]
                #print(j, prs[j-2])
                j=j+1
                
            else:
                if 'Inh' in prs[j-2][0] or 'Act' in prs[j-2][0]:
                    nor[index_nodes[prs[j-2][1]]][i] = nor[index_nodes[prs[j-2][1]]][i]*params[i][j]
                j=j+1

    ctr = 0
    ss_val = [[] for i in range(nodes)]   
    """Normalization with g,I1,I2,..K"""
    for i in range(len(ss)):
        if type(ss[i][0])==list or type(ss[i][0])==np.ndarray:
            for j in range(len(ss[i])):
                for k in range(int(ss[i][j][1])):
                    for m in range(nodes):
                        ss[i][j][nodes*k + m+2] = np.log2(((2**ss[i][j][nodes*k + m+2])/nor[m][int(ss[i][j][0]-1)]))
                        ss_val[m].append(ss[i][j][nodes*k + m +2])
                    ctr = ctr + 1
            
        else:
            for k in range(int(ss[i][1])):
                for m in range(nodes):
                    ss[i][nodes*k + m+2] = np.log2(((2**ss[i][nodes*k + m+2])/nor[m][int(ss[i][0]-1)]))
                    ss_val[m].append(ss[i][nodes*k + m +2])
                ctr = ctr + 1
            
                #print(n1[int(ss[i][j][0]-1)])

    for i in range(nodes):
        ss_val[i] = np.array(ss_val[i])
    avg = [ss_val[i].mean() for i in range(nodes)]
    std = [ss_val[i].std() for i in range(nodes)]
    zs = [zscore(ss_val[i]) for i in range(nodes)]
    
    return zs


name = 'Racipe_files1/2_MI' # Specify file location
name_nodes, core_nodes, sig_nodes, out_nodes, index_nodes = get_core_network(name) 
zscores = ss_zscore(name, core_nodes, index_nodes)
zscores = np.array(zscores)
"""
For fig 2, we plotted the heatmaps of zscores. 
"""

# Prints correlation values
print(np.corrcoef(zscores))
"""
We used the the scipy.stats.spearmanr function to obtain the spearman correlation 
for figure 4, 5 and 6. ANOVA analysis was done to obtain the error bars.
"""