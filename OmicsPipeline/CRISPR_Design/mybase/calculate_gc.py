def calculate_gc_content(dna_sequence):
    """
    计算DNA序列的GC含量。
    返回:
    GC含量（百分比）。
    """
    dna_sequence = dna_sequence[:-3]  ## 去除NGG序列
    count_g = dna_sequence.count('G')
    count_c = dna_sequence.count('C')
    total_nucleotides = len(dna_sequence)
    gc_content = round((count_g + count_c) / total_nucleotides * 100,2)
    return gc_content
