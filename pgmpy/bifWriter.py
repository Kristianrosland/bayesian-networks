import numpy as np
from itertools import product

def write_model(output_file, model):
    f = open(output_file, "w+")
    f.write("network " + (model.name or "unknown") + " {\n}\n")

    # Write variables 
    '''
        variable asia {
            type discrete [ 2 ] { yes, no };
        }
    '''
    for node in model.nodes(): 
        c = model.get_cardinality(node)
        possible_values = ", ".join(map(str, range(c)))
        f.write("variable " + node + " {\n")
        f.write("\ttype discrete [ {} ] ".format(c))
        f.write("{ " + possible_values + " };\n")
        f.write("}\n")


    # Write probabilities (CPDs)
    '''
    probability ( tub | asia ) {
        (yes) 0.05, 0.95;
        (no) 0.01, 0.99;
    }
    '''
    for cpd in model.get_cpds():
        values = cpd.get_values()
        f.write("probability ( {} ".format(cpd.variable))
        if len(cpd.variables) == 1:
            flattened_values = [v for val in values for v in val] # [[1],[2]] => [1,2]
            table_values = ", ".join(map(str, flattened_values))
            f.write(") {\n")
            f.write("\ttable {}\n".format(table_values))
            f.write("}\n")
        else:
            parents = ", ".join(cpd.variables[1:])
            f.write("| " + parents + " ) {\n")
            f.write(get_conditional_table(cpd))
            f.write("}\n")

    f.close()

def get_conditional_table(cpd):
    result = ""
    cardinalities = cpd.cardinality[1:]

    cardinality_ranges = map(lambda c: range(c), cardinalities)
    probabilities = np.dstack(cpd.get_values())[0]
    params = product(*cardinality_ranges)

    for idx, param in enumerate(params):
        p = probabilities[idx]
        if (len(param) == 1):
            param = "({})".format(param[0])
        result += "\t{} {};\n".format(param, ", ".join(map(str, p)))
    
    return result
