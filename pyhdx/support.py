def get_reduced_blocks(k_series, max_combine=2, max_join=5):
    block_length = list(k_series.cov.block_length.copy())

    i = 0
    subblock = 0
    while i < len(block_length):
        curr_size = block_length[i]
        if curr_size <= max_combine:
            del block_length[i]
            subblock += curr_size
        elif subblock != 0:
            # if block_length > 1:
            block_length.insert(i, subblock)
            subblock = 0
            i += 1
        else:
            i += 1
    if subblock != 0:
        block_length.append(subblock)

    changes = True
    while changes:
        i = 0
        changes = False
        while i < len(block_length):
            curr_size = block_length[i]
            if curr_size < max_join:
                changes = True
                if i == 0:  # beginning of list
                    block_length[i + 1] += curr_size
                    del block_length[i]
                elif i == len(block_length) - 1:  # end of the list
                    block_length[i - 1] += curr_size
                    del block_length[i]
                else:
                    if block_length[i - 1] < block_length[i + 1]:
                        block_length[i - 1] += curr_size
                    else:
                        block_length[i + 1] += curr_size
                    del block_length[i]
            else:
                i += 1

    return block_length

def get_constant_blocks(k_series, block_size=10, initial_block=5):
    num_repeats = (k_series.cov.prot_len - initial_block) // block_size
    remainder = (k_series.cov.prot_len - initial_block) % block_size

    blocks = [initial_block] + [block_size] * num_repeats
    if remainder:
        blocks += [remainder]

    return blocks


def reduce_inter(args):
    """
    #  https://github.com/brentp/interlap/blob/3c4a5923c97a5d9a11571e0c9ea5bb7ea4e784ee/interlap.py#L224
    # MIT Liscence
    >>> reduce_inter([(2, 4), (4, 9)])
    [(2, 4), (4, 9)]
    >>> reduce_inter([(2, 6), (4, 10)])
    [(2, 10)]
    """
    if len(args) < 2: return args
    args.sort()
    ret = [args[0]]
    for next_i, (s, e) in enumerate(args, start=1):
        if next_i == len(args):
            ret[-1] = ret[-1][0], max(ret[-1][1], e)
            break

        ns, ne = args[next_i]
        if e > ns or ret[-1][1] > ns:
            ret[-1] = ret[-1][0], max(e, ne, ret[-1][1])
        else:
            ret.append((ns, ne))
    return ret
