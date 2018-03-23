def add_suffix(seq_records):
    suffix = "AAAAAAAAAAGGGGGGGGGG"
    out = []
    for r in seq_records:
        out.append(r + suffix)
    return out
