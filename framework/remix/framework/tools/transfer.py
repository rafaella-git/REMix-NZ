import networkx as nx
import numpy as np
import pandas as pd

idx = pd.IndexSlice


def build_cmat(inc):
    """Apply Kirchhoff's loop rule to alternating-current networks.

    Parameters
    ----------
    inc : array
        incidence matrix.

    Returns
    -------
    pandas.Dataframe
        Kirchhoff's loop matrix (cycle matrix), shows which links belong to
        which loops.
    """
    df = pd.DataFrame(inc)
    df.index = pd.MultiIndex.from_tuples(df[0]).swaplevel(i=0, j=1)
    df = df[1].reset_index().replace({-1: "start", 1: "end"})
    df = df.set_index(["level_2", "level_3", "level_0", 1]).unstack()
    df.columns = df.columns.get_level_values(1)

    C_mat = []
    for year in list(set(df.index.get_level_values(0))):
        for seg in list(set(df.index.get_level_values(1))):
            df_seg = df.loc[idx[year, seg, :], idx[:]]
            network = nx.from_pandas_edgelist(df_seg, source="start", target="end")
            network_d = nx.from_pandas_edgelist(
                df_seg, source="start", target="end", create_using=nx.DiGraph
            )
            cycles = list(
                nx.minimum_cycle_basis(network)
            )  # We want the least amount of edges

            for i, c in enumerate(cycles):
                c_o = sorted(nx.find_cycle(network_d.subgraph(c), orientation="ignore"))
                for j in c_o:
                    e = str(
                        df_seg[
                            np.logical_and.reduce(
                                [
                                    df_seg["start"].isin(j[0:2]),
                                    df_seg["end"].isin(j[0:2]),
                                ]
                            )
                        ].index.values[0][2]
                    )
                    if c_o[0][2] == "forward":
                        C_mat.append(
                            (
                                ("c{}".format(i + 1), e, year, seg),
                                +1 if j[2] == "forward" else -1,
                            )
                        )
                    else:
                        C_mat.append(
                            (
                                ("c{}".format(i + 1), e, year, seg),
                                +1 if j[2] != "forward" else -1,
                            )
                        )
    return C_mat
