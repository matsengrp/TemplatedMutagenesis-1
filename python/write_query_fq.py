from Bio.Seq import Seq


def write_query(queries, query_names, seed_starts, seed_lens, file_prefix):
    """Writes the set of queries that we use as input to bowtie2

    """
    # make the files
    f_fw = open(file_prefix + "_fw.fq", "w+")
    f_rc = open(file_prefix + "_rc.fq", "w+")
    for (query, query_name, seed_start, seed_len) in zip(queries, query_names, seed_starts, seed_lens):
        # make the queries
        query_fw = query[seed_start:]
        query_rc = str(Seq(query[:seed_start + seed_len]).reverse_complement())
        # write the query names
        formatted_query_name = "@{}\n".format(query_name)
        f_fw.write(formatted_query_name)
        f_rc.write(formatted_query_name)
        # write the query sequences
        f_fw.write(query_fw + "\n+\n")
        f_rc.write(query_rc + "\n+\n")
        # write the quality scores
        f_fw.write("~" * len(query_fw) + "\n")
        f_rc.write("~" * len(query_rc) + "\n")
    # close the files
    f_fw.close()
    f_rc.close()
