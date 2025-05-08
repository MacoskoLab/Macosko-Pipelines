import argparse
import os

import pandas as pd
from xml.etree import ElementTree
import yaml


def abs_path(root: str, path: str) -> str:
    """
    Get the absolute path of a file.  If the path is absolute, return it.  If the path is relative, return the
    absolute path by joining the root and path.
    """
    if path.startswith("/"):
        return path
    if path.startswith("http"):
        return path
    return (root if root.endswith("/") else root + "/") + path


def get_args():
    parser = argparse.ArgumentParser(description="Generate sample sheet to convert bcl files to fastq files")
    indexes_group = parser.add_mutually_exclusive_group(required=True)
    cycles_group = parser.add_mutually_exclusive_group(required=True)
    indexes_group.add_argument(
        "-b",
        "--barcodes",
        help="Tab separated file with a header of columns SampleName, IndexBarcode1, IndexBarcode2, and optional Lane",
    )
    indexes_group.add_argument(
        "-i",
        "--indexes",
        help="YAML/JSON file containing {lane:[index,...], lane:[index,...], ...}",
    )
    cycles_group.add_argument(
        "-x",
        "--run-info",
        help="The RunInfo.xml from the bcl directory containing the list of cycles.",
    )
    cycles_group.add_argument(
        "-c",
        "--cycles",
        help="Comma-separated list of cycles for each read in the form of int,int,int,int",
    )
    parser.add_argument(
        "-l",
        "--lane",
        type=int,
        default=0,
        help="Lane number to subset the sample sheet to. Set to 0 to include all lanes. Default: %(default)s",
    )
    parser.add_argument(
        "-r",
        "--indexes_root",
        default="https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/bcl-convert/indexes/",
        help="Path to the directory containing csv files.  Default: %(default)s",
    )
    parser.add_argument(
        "--nn",
        default="SI-NN.csv",
        help="Absolute or relative path to the SN-NN.csv.  Default: %(default)s",
    )
    parser.add_argument(
        "--nt",
        default="SI-NT.csv",
        help="Absolute or relative path to the SI-NT.csv.  Default: %(default)s",
    )
    parser.add_argument(
        "--tt",
        default="SI-TT.csv",
        help="Absolute or relative path to the SI-TT.csv.  Default: %(default)s",
    )
    parser.add_argument(
        "--nd",
        default="ND7.csv",
        help="Absolute or relative path to the ND7.csv.  Default: %(default)s",
    )
    parser.add_argument(
        "-s",
        "--sample-sheet",
        required=True,
        help="Path to output the sample sheet file",
    )

    args = parser.parse_args()

    if args.cycles is not None:
        cycles = args.cycles.split(",")
        if len(cycles) != 4:
            parser.error("--cycles must contain 4 comma-separated int values")
        try:
            args.cycles_list = [int(cycle) for cycle in cycles]
        except ValueError:
            parser.error("--cycles must contain 4 comma-separated int values")

    if args.lane < 0 or args.lane > 8:
        parser.error("--lane must be between 0 and 8")

    return args


def barcodes_to_df(barcodes: pd.DataFrame) -> pd.DataFrame:
    """
    Read the barcode indexes from a tsv returning a dataframe with columns Lane, Sample_ID, index, and index2.
    """
    assert "SampleName" in barcodes.columns
    assert "IndexBarcode1" in barcodes.columns
    assert "IndexBarcode2" in barcodes.columns
    barcodes.rename(
        columns={"SampleName": "Sample_ID", "IndexBarcode1": "index", "IndexBarcode2": "index2"},
        inplace=True,
    )
    # In 'SampleName' replace anything not alphanumeric, -, or _ with '_'
    barcodes["Sample_ID"] = barcodes["Sample_ID"].str.replace(r"[^a-zA-Z0-9-_]", "_", regex=True)
    if not "Lane" in barcodes.columns:
        # For lane in 1 to 8, create a copy of the barcodes for each lane
        barcodes = pd.concat([barcodes] * 8, ignore_index=True)
        barcodes["Lane"] = [i + 1 for i in range(8) for _ in range(len(barcodes) // 8)]
    return barcodes[["Lane", "Sample_ID", "index", "index2"]]


def indexes_to_df(
        indexes: dict[int, list[str]],
        nn: pd.DataFrame,
        nt: pd.DataFrame,
        tt: pd.DataFrame,
        nd: pd.DataFrame,
) -> pd.DataFrame:
    """
    Convert the indexes dictionary to a DataFrame using the SI-NN, SI-NT, SI-TT, and ND7 csv files.
    """
    rows = []
    assert type(indexes) == dict
    for lane in sorted(indexes.keys()):
        assert lane in [1, 2, 3, 4, 5, 6, 7, 8]
        assert type(indexes[lane]) == list and len(indexes[lane]) == len(set(indexes[lane]))
        for index in indexes[lane]:
            if index in nd["I7_Index_ID"].values:
                row = nd.loc[nd['I7_Index_ID'] == index]
                seq1 = row["index"].iloc[0]
                rows.append([lane, index, seq1 + "AT", "ACCGAGATCT"])
            elif index in tt["index_name"].values:
                row = tt.loc[tt['index_name'] == index]
                seq1 = row["index(i7)"].iloc[0]
                seq2 = row["index2_workflow_a(i5)"].iloc[0]
                rows.append([lane, index, seq1, seq2])
            elif index in nt["index_name"].values:
                row = nt.loc[nt['index_name'] == index]
                seq1 = row["index(i7)"].iloc[0]
                seq2 = row["index2_workflow_a(i5)"].iloc[0]
                rows.append([lane, index, seq1, seq2])
            elif index in nn["index_name"].values:
                row = nn.loc[nn['index_name'] == index]
                seq1 = row["index(i7)"].iloc[0]
                seq2 = row["index2_workflow_a(i5)"].iloc[0]
                rows.append([lane, index, seq1, seq2])
            else:
                raise IndexError(f"ERROR: index {index} not found")

    df = pd.DataFrame(rows, columns=['Lane', 'Sample_ID', 'index', 'index2'])

    # If whole flowcell is 8bp, remove the buffer padding
    if df['index2'].eq('ACCGAGATCT').all():
        df['index'] = df['index'].str[:8]
        df['index2'] = df['index2'].str[:8]

    assert len(df) == len(df.drop_duplicates())

    return df


def read_run_info_cycles(run_info: str) -> list[int]:
    """
    Read the RunInfo.xml file and return the cycles for each read.
    """
    tree = ElementTree.parse(run_info)
    root = tree.getroot()
    cycles = [int(read.get("NumCycles")) for read in root.findall('Run/Reads/Read')]
    return cycles


def read_indexes(indexes_path: str) -> dict[int, list[str]]:
    """
    Read the indexes file and return a dictionary of lane numbers and their corresponding indexes.
    """
    with open(indexes_path, "r") as f:
        indexes = yaml.safe_load(f)
    return indexes


def write_sample_sheet(
        samplesheet: pd.DataFrame,
        cycles: list[int],
        sheet_path: str,
) -> None:
    """
    Write the sample sheet to a file.
    """
    sheet_path_dir = os.path.dirname(sheet_path)
    if sheet_path_dir:
        os.makedirs(sheet_path_dir, exist_ok=True)
    with open(sheet_path, 'w') as f:
        f.write("[Settings]\n")
        f.write("CreateFastqForIndexReads,0\n")  # default: 0
        f.write("NoLaneSplitting,false\n")  # default: false
        i1len = len(samplesheet["index"][0])
        i2len = len(samplesheet["index2"][0])
        assert type(cycles) == list and all(type(c) == int for c in cycles) and len(cycles) == 4
        assert cycles[0] > 0 and cycles[1] >= i1len and cycles[2] >= i2len and cycles[3] > 0
        r1 = f"Y{cycles[0]}"
        r2 = f"I{i1len}" + ("" if cycles[1] == i1len else f"N{cycles[1] - i1len}")
        r3 = f"I{i2len}" + ("" if cycles[2] == i2len else f"N{cycles[2] - i2len}")
        r4 = f"Y{cycles[3]}"
        f.write(f"OverrideCycles,{r1};{r2};{r3};{r4}\n")
        f.write("\n")

        f.write("[Data]\n")
    samplesheet.to_csv(sheet_path, mode='a', index=False)


def run() -> None:
    args = get_args()

    if args.barcodes:
        print(f"Barcodes: {args.barcodes}")
    else:
        print(f"Indexes: {args.indexes}")
    if args.run_info:
        print(f"RunInfo: {args.run_info}")
    else:
        print(f"Cycles: {args.cycles_list}")
    print(f"Indexes root: {args.indexes_root}")
    print(f"NN: {args.nn}")
    print(f"NT: {args.nt}")
    print(f"TT: {args.tt}")
    print(f"ND: {args.nd}")
    print(f"Sample sheet: {args.sample_sheet}")

    # Get the cycles from the run info or the command line
    if args.run_info:
        cycles_list = read_run_info_cycles(args.run_info)
    else:
        cycles_list = args.cycles_list

    # Get the barcode indexes from barcodes tsv or the indexes yaml
    if args.barcodes:
        barcodes_df = pd.read_csv(args.barcodes, sep="\t")
        samplesheet_df = barcodes_to_df(barcodes_df)
    else:
        with open(args.indexes, "r") as f:
            indexes_dict = yaml.safe_load(f)
        nn = pd.read_csv(abs_path(args.indexes_root, args.nn))
        nt = pd.read_csv(abs_path(args.indexes_root, args.nt))
        tt = pd.read_csv(abs_path(args.indexes_root, args.tt))
        nd = pd.read_csv(abs_path(args.indexes_root, args.nd))
        samplesheet_df = indexes_to_df(indexes_dict, nn, nt, tt, nd)

    if args.lane > 0:
        samplesheet_df = samplesheet_df[samplesheet_df["Lane"] == args.lane].reset_index(drop=True)

    write_sample_sheet(samplesheet_df, cycles_list, args.sample_sheet)


if __name__ == "__main__":
    run()
