import pandas as pd
import argparse


def options():
    parser = argparse.ArgumentParser()
    general = parser.add_argument_group('General info about the data.')
    general.add_argument('db', type=str, default=None,
                         help='Path to pickle file to operate on.')
    general.add_argument('--exp', type=str, default=None,
                         help='Experiment name')
    general.add_argument('--dish', type=str, default=None,
                         help='Dish name')
    general.add_argument('--scan', type=str, default=None,
                         help='Scan name')
    general.add_argument('--colname', type=str, default=None,
                         help='Column name to add.')
    general.add_argument('--val', type=str, default=None,
                         help='Value that should go in column to add.')
    general.add_argument('--unique', action='store_true',
                         help='If set, will enforce that --val is unique in --colname.')
    return parser.parse_args()


def add_col(exp, dish, scan, colname, colval, data, unique=False):
    if unique:
        if colname in data.columns:
            colidx = data.columns.get_loc(colname)
            if colval in data.iloc[:, colidx].unique():
                raise ValueError(f'Value {colval} already exists in column {colname} in.')
    if data[(data.experiment == exp) & (data.dish == dish) & (data.scan == scan)].empty:
        raise ValueError(f'No entry for {exp}_{dish}_{scan}.')
    data.loc[(data.experiment == exp) &
             (data.dish == dish) &
             (data.scan == scan),
             colname] = colval
    return data


if __name__ == "__main__":
    args = options()
    data = pd.read_pickle(args.db)
    data = add_col(args.exp, args.dish, args.scan,
                   args.colname, args.val, data, args.unique)
    data.to_pickle(args.db)
