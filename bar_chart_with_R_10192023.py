import argparse
import csv
import matplotlib.pyplot as plt
import numpy as np
import textwrap
import matplotlib.colors as mcolors
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from rpy2.robjects.conversion import localconverter
import pandas as pd
import copy

def parse_arguments():
    parser = argparse.ArgumentParser(description='Bar chart generation')
    parser.add_argument('-f', '--filename', help='Path to the CSV file')
    parser.add_argument('-t', '--test', help='Default is LSD. If test flag is set to be turkey, turkey test will be conducted.')
    return parser.parse_args()

def read_csv_file(filename):
    with open(filename, "r") as csvfile:
        csvreader = csv.reader(csvfile)
        y_axis_name = next(csvreader)[1]
        group_names_list = next(csvreader)
        group_names_list_nonull = [x for x in group_names_list if x]
        if len(group_names_list_nonull) == 1:
            group_name = group_names_list_nonull[0]
            data = [row for row in csvreader]
        elif len(group_names_list_nonull) == 2:
            group_subgroup_data = {}
            prev_group_name = ""
            for row in csvreader:
                if row[0]:
                    group_name = row[0]
                    prev_group_name = row[0]
                    subgroup_name = row[1]
                    subgroup_data = [float(data) for data in row[2:] if data]
                else:
                    group_name = prev_group_name
                    subgroup_name = row[1]
                    subgroup_data = [float(data) for data in row[2:] if data]
                if group_name not in group_subgroup_data:
                    group_subgroup_data[group_name] = {}
                group_subgroup_data[group_name][subgroup_name] = subgroup_data
            data = group_subgroup_data
        return y_axis_name, data

def run_r_test(data, testname):
    # group_one_names = list(data.keys())[:2]
    r_data = copy.deepcopy(data)
    r_df = []
    dropped_key = []
    # Check if r_data is a dictionary of dictionaries
    if not isinstance(r_data, dict):
        for group, values in r_data:
            temp = values
            temp.insert(0, str(group))
            r_df.append(temp)

    else:
        for group, subgroups in r_data.items():
            for subgroup, values in subgroups.items():
                temp = values
                temp.insert(0, str(group) + ',' + str(subgroup))
                r_df.append(temp)
    
    data = r_df

    pandas2ri.activate()
    base = importr("base")
    stats = importr("stats")
    agricolae = importr("agricolae")

    data_long_list = [
        (row[0], i + 1, value) for row in data for i, value in enumerate(row[1:])
    ]

    data_long_df = base.data_frame(
        Group=ro.StrVector([row[0] for row in data_long_list]),
        Value=ro.FloatVector([row[2] for row in data_long_list]),
    )

    robjects.globalenv['data_long_df'] = data_long_df

    if testname == "turkey":
        r_code = '''
        model <- aov(Value ~ Group, data = data_long_df)
        cv.model(model)
        df <- df.residual(model)
        MSerror <- deviance(model) / df
        HSD.test(model, "Group", console = TRUE)
        '''
    else:
        r_code = '''
        model <- aov(Value ~ Group, data = data_long_df)
        cv.model(model)
        df <- df.residual(model)
        MSerror <- deviance(model) / df
        LSD.test(model, "Group",  p.adj="bonferroni", console = TRUE)
        '''

    result = robjects.r(r_code)
    with localconverter(ro.default_converter + pandas2ri.converter):
          pd_from_r_df = ro.conversion.rpy2py(result[4])

    return pd_from_r_df

def generate_bar_chart(y_axis_name, data, r_test_results):
    subgroups = np.unique([s for d in data.values() for s in d.keys()])
    num_groups = len(data)

    patterns = ['/', '\\', 'x', '*', '.', '+', 'o', '|', 'O', '-',]
    colors = list(mcolors.TABLEAU_COLORS.values())

    means = np.zeros((len(subgroups), num_groups))
    stds = np.zeros((len(subgroups), num_groups))
    for i, subgroup in enumerate(subgroups):
        for j, (group, values) in enumerate(data.items()):
            if subgroup in values:
                means[i, j] = np.mean(values[subgroup])
                stds[i, j] = np.std(values[subgroup], ddof=1)

    fig, ax = plt.subplots(figsize=(10, 5))
    width = 1 / (num_groups + 1)
    for j, group in enumerate(data):
        ax.bar(np.arange(len(subgroups)) + j * width, means[:, j], width=width, yerr=stds[:, j], label=group, alpha=0.75, capsize=10, 
               color=colors[j], hatch=patterns[j], ecolor='black', error_kw=dict(lw=1, capthick=1, capsize=5))

    for i, subgroup in enumerate(subgroups):
        for j, (group, values) in enumerate(data.items()):
            if subgroup in values:
                key = str(group) +','+ str(subgroup)
                label = r_test_results.loc[key]['groups']
                value = means[i, j]
                if value != 0:
                    ax.text(i + j * width, value + stds[i, j], label, ha='center', va='bottom')

    ax.set_xticks(np.arange(len(subgroups)) + (num_groups - 1) * width / 2)

    subgroups_wrapped = [textwrap.fill(s, 10) for s in subgroups]
    ax.set_xticklabels(subgroups_wrapped)
    ax.set_ylabel(y_axis_name)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=2)
    ax.annotate('Fast Bar Chart', (2.5, 10), alpha=0.2, ha="center", va="center", fontsize=30,  rotation=30, color="gray")
    manager = plt.get_current_fig_manager()
    manager.full_screen_toggle()
    plt.tight_layout()
    plt.show()

def main():
    args = parse_arguments()
    y_axis_name, data = read_csv_file(args.filename)
    r_test_results = run_r_test(data, args.test)
    generate_bar_chart(y_axis_name, data, r_test_results)

if __name__ == '__main__':
    main()