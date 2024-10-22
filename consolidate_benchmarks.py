#!/usr/bin/env python3

import os
import pandas as pd
import sys
import xlsxwriter
from xlsxwriter.utility import xl_rowcol_to_cell

sort_priority = {
    "640_AES": 0,
    "640_SHAKE": 1,
    "976_AES": 2,
    "976_SHAKE": 3,
    "1344_AES": 4,
    "1344_SHAKE": 5,
    "lightsaber": 6,
    "saber": 7,
    "firesaber": 8,
    "stack": 9,
    "mmap": 10,
    "ref": 11,
    "opt": 12,
    "neon": 13,
    "opt_amx": 14,
    "amx_polymul": 15,
    "amx_matmul": 16,
}


def my_sort(x):
    return sort_priority.get(x, 1337)


def my_sort_pd(series):
    return series.apply(lambda x: sort_priority.get(x, 1337))


def write_benchmarks(scheme, worksheet, df, columns):
    center_h = workbook.add_format({"align": "center"})
    center_h_bold = workbook.add_format({"align": "center", "bold": True})
    center_h_bold_2places = workbook.add_format(
        {"align": "center", "bold": True, "num_format": "0.00"}
    )
    center_hv = workbook.add_format({"align": "center", "valign": "vcenter"})

    # "Scheme or operation",
    # "Parameter set",
    # "Memory allocation",
    # "Implementation",
    # "Variant",
    # "Operation",
    # "Cycle count",

    match scheme:
        case "frodokem":
            benchmark_names = [
                ("Key generation", "crypto_kem_keypair"),
                ("Encapsulation", "crypto_kem_enc"),
                ("Decapsulation", "crypto_kem_dec"),
                ("Encapsulation 4x", "crypto_kem_enc4x"),
                ("Decapsulation 4x", "crypto_kem_dec4x"),
                ("A*s + e", "FrodoKEM A*s + e"),
                ("A*s + e matmul", "FrodoKEM A*s + e (matmul only)"),
                ("s*A + e", "FrodoKEM s*A + e"),
                ("s*A + e matmul", "FrodoKEM s*A + e (matmul only)"),
                ("s*A + e 4x", "FrodoKEM s*A + e 4x"),
                ("s*A + e matmul 4x", "FrodoKEM s*A + e 4x (matmul only)"),
            ]
        case "saber":
            benchmark_names = [
                ("Key generation", "crypto_kem_keypair"),
                ("Encapsulation", "crypto_kem_enc"),
                ("Decapsulation", "crypto_kem_dec"),
                ("MatrixVectorMulRound", "MatrixVectorMulRound"),
            ]

    impls = {}

    for pset in [
        "640_AES",
        "640_SHAKE",
        "976_AES",
        "976_SHAKE",
        "1344_AES",
        "1344_SHAKE",
    ]:
        impls[pset] = {
            "ref": {
                "ref": "Reference",
                "opt": "Optimized",
                "neon": "Ours (NEON)",
                "opt_amx": "Ours (AMX optimized)",
            }
        }

    for pset in ["lightsaber", "saber", "firesaber"]:
        impls[pset] = {
            "BHK": {
                "neon": "[BHK21]",
                "amx_polymul": "Ours (polynomial multiplication)",
                "amx_matmul": "Ours (matrix multiplication)",
            }
        }

    for i in range(3):
        worksheet.merge_range(0, i, 1, i, columns[i + 1], center_hv)

    worksheet.merge_range(0, 3, 0, 3 + len(benchmark_names) - 1, "Operation", center_h)

    for i, names in enumerate(benchmark_names):
        worksheet.write(1, 3 + i, names[0], center_h)

    current_row = 2
    for pset in sorted(df["Parameter set"].unique(), key=my_sort):
        merge_pset_start = current_row

        df2 = df[df["Parameter set"] == pset]

        for memalloc in sorted(df2["Memory allocation"].unique(), key=my_sort):
            merge_memalloc_start = current_row

            df3 = df2[df2["Memory allocation"] == memalloc]

            for work in sorted(df3["Implementation"].unique(), key=my_sort):
                df4 = df3[df3["Implementation"] == work]

                for variant in sorted(df4["Variant"].unique(), key=my_sort):
                    impl = impls[pset][work][variant]
                    if not impl:
                        continue

                    worksheet.write(current_row, 2, impl, center_h)

                    df5 = df4[df4["Variant"] == variant]

                    for name, op in benchmark_names:
                        try:
                            cycle_count = int(
                                df5[df5["Operation"] == op]["Cycle count"].iloc[0]
                            )

                            worksheet.write(
                                current_row,
                                3 + benchmark_names.index((name, op)),
                                cycle_count,
                                center_h,
                            )
                        except IndexError:
                            pass

                    current_row += 1

            worksheet.merge_range(
                merge_memalloc_start, 1, current_row - 1, 1, memalloc, center_hv
            )

        worksheet.merge_range(current_row, 1, current_row, 2, "Speedup", center_h_bold)

        for j in range(3, 3 + len(benchmark_names)):
            formula = (
                "{=ROUND(MIN(IF(RIGHT("
                f"{xl_rowcol_to_cell(merge_pset_start, 2)}:"
                f'{xl_rowcol_to_cell(current_row - 1, 2)})<>")",'
                f"{xl_rowcol_to_cell(merge_pset_start, j)}:"
                f'{xl_rowcol_to_cell(current_row - 1, j)},""'
                "))/MIN(IF(RIGHT("
                f"{xl_rowcol_to_cell(merge_pset_start, 2)}:"
                f'{xl_rowcol_to_cell(current_row - 1, 2)})=")",'
                f"{xl_rowcol_to_cell(merge_pset_start, j)}:"
                f"{xl_rowcol_to_cell(current_row - 1, j)},"
                '"")),2)}'
            )

            worksheet.write(current_row, j, formula, center_h_bold_2places)

        worksheet.merge_range(merge_pset_start, 0, current_row, 0, pset, center_hv)

        current_row += 1

    worksheet.autofit()


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("Usage: python consolidate_benchmarks.py <CPU name>")
        sys.exit(1)

    path = os.path.join(os.getcwd(), f"speed_results_{sys.argv[1]}")

    if not os.path.isdir(path):
        print(f"Directory {path} does not exist.")
        sys.exit(1)

    all_files = [
        f
        for f in os.listdir(path)
        if os.path.isfile(os.path.join(path, f)) and ":" in f and f.endswith(".txt")
    ]

    data = []

    for file in all_files:
        row = file.split(".txt")[0].split(":")

        with open(os.path.join(path, file), "r", encoding="utf-8") as f:
            for line in f.read().splitlines():
                if line:
                    data.append(row + line.split(":"))

    columns = [
        "Scheme or operation",
        "Parameter set",
        "Memory allocation",
        "Implementation",
        "Variant",
        "Operation",
        "Cycle count",
    ]

    df = pd.DataFrame(data, columns=columns).sort_values(
        by=columns[:-2], key=my_sort_pd
    )

    matcher = {
        "frodokem": "(?:frodokem|matmul)",
        "saber": "(?:saber|matrixvectormulround)",
    }

    workbook = xlsxwriter.Workbook(
        os.path.join(f"speed_results_{sys.argv[1]}", "benchmarks.xlsx")
    )
    worksheet = {}
    for scheme in ["frodokem", "saber"]:  # ["frodokem", "frodokem4x", "saber"]:
        worksheet[scheme] = workbook.add_worksheet(scheme)

        write_benchmarks(
            scheme,
            worksheet[scheme],
            df[df["Scheme or operation"].str.contains(matcher[scheme])],
            columns,
        )

    workbook.close()
