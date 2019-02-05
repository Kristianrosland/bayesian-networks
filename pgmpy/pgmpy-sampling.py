import pandas as pd # pylint: disable=E0401
from pgmpy.sampling import BayesianModelSampling
from pgmpy.readwrite import BIFReader
from bifWriter import write_model

dataset_name = "alarm"

## FROM THE .bif FILE
reader = BIFReader("../data/true_{}.bif".format(dataset_name))
model = reader.get_model()
model.check_model()

sample_size = 100000
sampler = BayesianModelSampling(model)
sampled_dataframe = sampler.forward_sample(size=sample_size)

num_columns = len(sampled_dataframe.columns)
columns_translate_dict = {}
cardinalities = []

for idx in range(num_columns):
    column_name = sampled_dataframe.columns[idx]
    cardinalities.append(model.get_cardinality(column_name))
    columns_translate_dict[idx] = column_name

sampled_dataframe.columns = range(num_columns)

generated_data = open("../data/{}_{}.data".format(dataset_name, sample_size), "w+")
generated_data.write("{}\n".format(model.number_of_nodes()))
generated_data.write(" ".join(map(str, cardinalities)))
generated_data.write("\n")
generated_data.write(str(sample_size) + "\n")

for (_,row) in sampled_dataframe.iterrows():
    values = map(lambda col: row[col], sampled_dataframe.columns)
    generated_data.write("{}\n".format(" ".join(map(str, values))))

translate_file = open("../{}_translate.txt".format(dataset_name), "w+")
for x in range(num_columns):
    translate_file.write("{}: {}\n".format(x, columns_translate_dict[x]))

generated_data.close()
translate_file.close()

'''
# network_means = sampled_dataframe.mean().sort_index()
# network_means = network_means.map(lambda x: round(x, 2))
# print network_means

## FROM THE .data FILE
input_file = open("../data/alarm_100.data")
num_vars = int(input_file.readline())
num_values = map(int, input_file.readline().split(' '))
num_datapoints = int(input_file.readline())

datapoints = []
for line in input_file:
    datapoints.append(map(int, line.split(' ')))

dataframe = pd.DataFrame(datapoints)
actual_means = dataframe.mean()

## Log output
# log_file = open("../logs/pgmpy-sampling.log")
diffs = map(lambda (a,b): abs(round(a-b, 2)), zip(actual_means, network_means))
compared_means = pd.DataFrame(
    {
        'sampled': network_means.values,
        'actual': actual_means.values,
        'diff': diffs
    })

# print compared_means[['sampled', 'actual', 'diff']]
'''