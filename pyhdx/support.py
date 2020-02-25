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
