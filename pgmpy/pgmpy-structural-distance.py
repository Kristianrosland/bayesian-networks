import random
import argparse
import networkx as nx
from pgmpy.readwrite import BIFReader
from pgmpy.models import BayesianModel
    
parser = argparse.ArgumentParser()
parser.add_argument("-f", help="learnt network .bif file")
parser.add_argument("-t", help="true network .bif file")
parser.add_argument("-tr", help="translation file e.g. (0: HVPR)")
args = parser.parse_args()


def has_unordered_edge(node, edges):
    for (_, v, idx) in edges:
        if node == v and idx == -1:
            return True
    return False

'''
    Order according to Chickerings algorithm part 1 (1995)
'''
def order_edges(network):
    nodes = nx.topological_sort(network)
    reversed_nodes = nodes[::-1]
    edges = map(lambda (u,v): (u,v,-1), network.edges())
    index = 0

    while (index < len(edges)):
        # lowest ordered node that has an unordered edge incident (pointing?) to it
        y = next(node for node in nodes if has_unordered_edge(node, edges)) 

        edges_to_y = filter(lambda (u,v, idx): idx == -1 and v == y, edges)
        nodes_connected_to_y = map(lambda (u, v, idx): u, edges_to_y)
        
        # the highest ordered node x for which x -> y is not ordered
        x = next(node for node in reversed_nodes if node in nodes_connected_to_y)

        # update edge (u,v,-1) to (x, y, index)
        for idx in range(len(edges)):
            (u,v,i) = edges[idx]
            if u == x and v == y and i == -1:
                edges[idx] = (x,y,index)
                index += 1
        

    sorted_edges = sorted(edges, key=lambda (u,v,idx): idx)
    return map(lambda (u,v,_): (u,v,'unknown'), sorted_edges)


'''
    Labels edges either 'compelled' or 'reversible'
'''
def label_edges(network):
    ordered_edges = order_edges(network)

    while (len(filter(lambda (u, v, label): label == 'unknown', ordered_edges)) != 0):
        fixed_edges = []
        iteration_done = False

        (x, y, _) = next(edge for edge in ordered_edges if edge[2] == 'unknown')
        unknown_edges_into_y = filter(lambda(_, y_, label): y_ == y and label == 'unknown', ordered_edges)
        compelled_edges_to_x = filter(lambda(u,v,label): v == x and label == 'compelled', ordered_edges)

        for (w, _, _) in compelled_edges_to_x: # For every edge w -> x labeled 'compelled'
            if not w in network.get_parents(y): # If w is not a parent of y
                fixed_edges.append((x,y,'compelled')) #  Label x -> y with 'compelled'
                for (x_, y_, _) in ordered_edges: # and label every edge incident into y with 'compelled'
                    if y_ == y:
                        fixed_edges.append((x_, y_, 'compelled'))
                        iteration_done = True
            else:
                fixed_edges.append((w, y, 'compelled'))

        
        z_to_y = filter(lambda (z, y_, l_): z != x and y_ == y and not z in network.get_parents(x), ordered_edges)
        if not iteration_done:
            if len(z_to_y) != 0:  # if there exists an edge z -> y such that z != x and z is not a parent of x
                label = 'compelled'
            else:
                label = 'reversible'
            
            fixed_edges.append((x,y,label))
            for (x_, y_, _) in unknown_edges_into_y:
                fixed_edges.append((x, y, label))
        
        # update ordered edges
        for (u, v, l) in fixed_edges:
            for i in range(len(ordered_edges)):
                (u_, v_, _) = ordered_edges[i]
                if u == u_ and v == v_:
                    ordered_edges[i] = (u, v, l)

    return ordered_edges

def missing_edges(true_network, cpdag):
    count = 0
    for (u, v) in true_network.edges():
        if not any((x == u and y == v) or (x == v and y == u and label == 'reversible') for (x, y, label) in cpdag):
            count += 1

    return count

def added_edges(true_network, cpdag):
    count = 0
    for (x, y, label) in cpdag:
        if label == 'compelled':
            if not any(x == u and y == v for (u,v) in true_network.edges()):
                count += 1
        elif label == 'reversible':
            if not any((x == u and y == v) or (x == v and y == u) for (u,v) in true_network.edges()):
                count += 1
    return count

def incorrect_direction(true_network, cpdag):
    count = 0
    for (x, y, label) in cpdag:
        if label == 'compelled':
            if any(x == v and y == u for (u, v) in true_network.edges()):
                count += 1
    return count


reader = BIFReader(args.f)
learnt_network = reader.get_model()
learnt_network.check_model()

reader = BIFReader(args.t)
true_network = reader.get_model()
true_network.check_model()

# CPDAG / pattern of learnt network (equivalence class)
cpdag = label_edges(learnt_network)
if not args.tr is None:
    translation_file = open(args.tr, "r")
    translation_dict = {} 
    for line in translation_file:
        [number, variable] = line.strip('\n').split(": ")
        translation_dict[number] = variable
    cpdag = map(lambda (u,v,label): (translation_dict[u], translation_dict[v], label), cpdag)
    translated_edges = map(lambda (a,b): (translation_dict[a], translation_dict[b]) ,learnt_network.edges())


total_edges_true = true_network.number_of_edges()
total_edges_learnt = len(cpdag)
missing_edges = missing_edges(true_network, cpdag)
added_edges = added_edges(true_network, cpdag)
incorrect_edges = incorrect_direction(true_network, cpdag)
hamming_distance = missing_edges + added_edges + incorrect_edges


print "Total edges in learnt network: {}".format(total_edges_learnt)
print "Total edges in true network: {}\nMissing edges: {}\nAdded edges: {}\nIncorrect dir: {}\n".format(total_edges_true, missing_edges, added_edges, incorrect_edges)
print "Structural Hamming Distance: {}".format(hamming_distance)




'''
random_network = BayesianModel()
random_edges = []
variables = true_network.nodes()

for i in range(learnt_network.number_of_edges()):
    x = variables[random.randint(0, len(variables)-1)]
    y = variables[random.randint(0, len(variables)-1)]
    if x == y:
        i -= 1
    else:
        random_edges.append((x, y))

# random_network.add_edges_from(random_edges)

# network.add_edges_from([('A','C'),('B','D'),('C','E'),('D','E'),('E','F'),('E','G'),('H','G'), ('B', 'H')])
# network.add_edges_from([('diff', 'grade'), ('intel', 'grade'),('intel', 'SAT'), ('grade', 'letter'), ('SAT', 'letter'), ('SAT', 'spouse'), ('grade', 'spouse')])
# network.add_edges_from([('E', 'H'),('A','O'),('O','D'),('A','H'),('B','C'),('D','E'),('G','F'),('D','F'),('C','D'),('B','O')])
'''