import numpy as np
from pyhdx import Coverage, PeptideMasterTable
from numpy.lib.recfunctions import stack_arrays, append_fields


def gen_coverage(coverage, prot_length, pep_length):
    """

    :param coverage: target average number of peptides per residue
    :param length:
    :return:
    """

    num = int(prot_length*coverage / pep_length[0])
    peptides = []
    while len(peptides) < num:
        size = max(int(np.random.normal(*pep_length)), 3)
        print(size)
        middle = np.random.randint(0, prot_length+1)
        start = max(int(middle - size / 2), 0)
        end = min(int(middle + size / 2), prot_length)
        if (start, end) in peptides:
            print('present')
        else:
            peptides.append((start, end))


    peptides.sort()
    return peptides

def trace(length, section_len, corr, min_len=5):
    """generate uptake trace with given length consisting of sections with
    normally distbuted section_len (mean, sigma) who make jumps of corr (mean, sigma)

    """
    output = []

    initial = 100 * np.random.random()
    num = max(min_len, int(np.round(np.random.normal(*section_len))))
    output += [initial] * num
    while len(output) < length:
        dst = np.random.normal(*corr)
        sign = 2 * np.random.randint(2) - 1
        new = output[-1] + (dst * sign)
        if new > 100:
            new = 200 - new
        elif new < 0:
            new = -new

        num = max(min_len, int(np.round(np.random.normal(*section_len))))
        output += [new] * num

    return np.array(output[:length])


def to_rates(trace, kmin, kmax):
    """covert a trace of normalized values into kinetic rates in log space which are wihtin kmin and kmax"""
    lkmin = np.log10(kmin)
    lkmax = np.log10(kmax)

    x = (lkmax - lkmin) / 100
    log_rates = trace*x + lkmin
    rates = 10**log_rates

    return rates


def generate_data(peptides, sequence, timepoints, rates, state='state1'):
    """
    Generate HDX data array

    peptides: list of (start, stop)
    sequence: string of characters
    timepoints: list of timepoints, time units the same as rates
    rates: GT rates of exchange for peptides
    """

    start, end = np.array(peptides).T
    size = np.max(end-start) + 1

    # Generate a data array to create a coverage object
    dtype = [('start', int), ('end', int), ('exposure', float), ('state', f'U{len(state)}'), ('sequence', f'U{size}')]
    data = np.empty(len(start), dtype=dtype)
    data['start'] = start
    data['end'] = end
    data['exposure'][:] = 0.1
    data['state'][:] = 'state1'
    data['sequence'] = list([sequence[s-1:e] for s, e in zip(start, end)])

    pmt = PeptideMasterTable(data, drop_first=0, ignore_prolines=False, remove_nan=False)
    cov = Coverage(pmt.data)
    #crop rates to the right size
    prot_rates = rates[:cov.prot_len]

    arrays = []
    for t in timepoints:
        uptake = 1 - np.exp(-prot_rates * t)
        uptake *= 100  # pertange uptake per residue for time t
        scores = cov.X.dot(uptake)  # scores (%)

        full_dtype = dtype + [('scores', float)]
        full_data = np.empty(len(start), dtype=full_dtype)
        full_data['start'] = start
        full_data['end'] = end
        full_data['exposure'][:] = t
        full_data['state'][:] = 'state1'
        full_data['sequence'] = data['sequence']
        full_data['scores'] = scores

        arrays.append(full_data)

    stacked = stack_arrays(arrays, usemask=False, autoconvert=True)

    return cov, stacked, prot_rates
