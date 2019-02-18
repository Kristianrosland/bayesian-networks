import pandas as pd # pylint: disable=E0401

def get_translate_dict(translate_file):
    f = open(translate_file, "r")
    translation_dict = {} 
    for line in f:
        line = line.replace(": ", ":")
        [number, variable] = line.strip('\n').split(":")
        translation_dict[number] = variable
    return translation_dict

def import_and_translate_data(input_file, translate_file):
    f = open(input_file, "r")

    variables = f.readline().split(' ')
    translate_dict = get_translate_dict(translate_file)
    columns = map(lambda x: translate_dict[x], variables)

    data_points = []

    for line in f:
        data_points.append(map(int, line.split(' ')))

    return pd.DataFrame(data_points, columns=columns)
    