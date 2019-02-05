import networkx as nx
from pgmpy.readwrite import BIFReader

## FROM THE .bif FILE
reader = BIFReader("../data/true_alarm.bif")
model = reader.get_model()
model.check_model()

G = nx.DiGraph()
G.add_nodes_from(model.nodes())
G.add_edges_from(model.edges())

nx.write_gml(G, "../testing.gml")