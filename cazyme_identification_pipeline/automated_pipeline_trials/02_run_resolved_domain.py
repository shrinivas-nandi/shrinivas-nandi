#02_run_resolve_domain.py
################ Compare and look for overlapping domains on the same proteins ######################
import pandas as pd
import sys

# Column positions (0-indexed)
COL_QUERY    = 3   # protein name
COL_IEVALUE  = 12  # i_evalue
COL_ALI_FROM = 17  # ali_from
COL_ALI_TO   = 18  # ali_to
COL_COVERAGE = 22  # domain_coverage


def calc_overlap(ali_from1, ali_to1, ali_from2, ali_to2):
    overlap_start = max(ali_from1, ali_from2)
    overlap_end   = min(ali_to1, ali_to2)
    overlap_len   = max(0, overlap_end - overlap_start + 1)

    shorter = min(ali_to1 - ali_from1 + 1, ali_to2 - ali_from2 + 1)
    if shorter == 0:
        return 0
    return overlap_len / shorter


def resolve_overlaps(df, overlap_threshold=0.30):
    kept = []

    for protein, group in df.groupby(COL_QUERY):
        group = group.reset_index(drop=True)

        if len(group) == 1:
            kept.append(group)
            continue

        to_drop = set()
        rows = list(group.iterrows())

        for i, (idx1, row1) in enumerate(rows):
            if idx1 in to_drop:
                continue

            for j, (idx2, row2) in enumerate(rows):
                if i >= j:
                    continue
                if idx2 in to_drop:
                    continue

                overlap = calc_overlap(
                    row1[COL_ALI_FROM], row1[COL_ALI_TO],
                    row2[COL_ALI_FROM], row2[COL_ALI_TO]
                )

                if overlap >= overlap_threshold:
                    # lower i_evalue wins
                    if row1[COL_IEVALUE] < row2[COL_IEVALUE]:
                        to_drop.add(idx2)
                    elif row2[COL_IEVALUE] < row1[COL_IEVALUE]:
                        to_drop.add(idx1)
                    else:
                        # tie — higher coverage wins
                        if row1[COL_COVERAGE] >= row2[COL_COVERAGE]:
                            to_drop.add(idx2)
                        else:
                            to_drop.add(idx1)

        kept.append(group[~group.index.isin(to_drop)])

    return pd.concat(kept).reset_index(drop=True)


def main(input_file):
    df = pd.read_csv(input_file, sep="\t", header=None)

    # convert numeric columns
    for col in [COL_IEVALUE, COL_ALI_FROM, COL_ALI_TO, COL_COVERAGE]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    print(f"Total rows before: {len(df)}")

    df_resolved = resolve_overlaps(df, overlap_threshold=0.30)

    print(f"Total rows after:  {len(df_resolved)}")
    print(f"Rows removed:      {len(df) - len(df_resolved)}")

    # automatically generate output filename
    output_file = input_file.replace(".tsv", "_coverage_resolved_domain_output.tsv")

    df_resolved.to_csv(output_file, sep="\t", index=False, header=False)

    print(f"Output written to: {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python resolve_overlaps.py input.tsv")
        sys.exit(1)

    main(sys.argv[1])
