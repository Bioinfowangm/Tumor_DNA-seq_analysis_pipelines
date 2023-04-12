"""bcbio_freebayes_somatic_filter.py
Usage: bcbio_freebayes_somatic_filter.py [inVcf] [outVcf]
Will accept 'stdin' for inVcf.
Will accept 'stdout' for outVcf.
Extracts somatic variants from FreeBayes output, using bcbio filters.
Modified from https://github.com/chapmanb/bcbio-nextgen/blob/4fe770cc1343f8e1a3f3fab1771bad13eb94df7a/bcbio/variation/freebayes.py
"""

import os
import sys

if len(sys.argv) == 3:
    if sys.argv[1] == "stdin":
        inVcf = sys.stdin
    else:
        inVcf = open(sys.argv[1], 'r')
    if sys.argv[2] == "stdout":
        outVcf = sys.stdout
    else:
        outVcf = open(sys.argv[2], 'w')
else:
    print >> sys.stderr, \
        __doc__,
    sys.exit (1)

def _check_lods(parts, tumor_thresh, normal_thresh):
    """Ensure likelihoods for tumor and normal pass thresholds.

    Skipped if no FreeBayes GL annotations available.
    """
    try:
        gl_index = parts[8].split(":").index("GL")
    except ValueError:
        return True
    try:
        tumor_gls = [float(x) for x in parts[9].split(":")[gl_index].split(",")]
        tumor_lod = max(tumor_gls[i] - tumor_gls[0] for i in range(1, len(tumor_gls)))
    # No GL information, no tumor call (so fail it)
    except (IndexError, ValueError):
        tumor_lod = -1.0
    try:
        normal_gls = [float(x) for x in parts[10].split(":")[gl_index].split(",")]
        normal_lod = min(normal_gls[0] - normal_gls[i] for i in range(1, len(normal_gls)))
    # No GL inofmration, no normal call (so pass it)
    except (IndexError, ValueError):
        normal_lod = normal_thresh
    return normal_lod >= normal_thresh and tumor_lod >= tumor_thresh

def _check_freqs(parts):
    """Ensure frequency of tumor to normal passes a reasonable threshold.

    Avoids calling low frequency tumors also present at low frequency in normals,
    which indicates a contamination or persistent error.
    """
    thresh_ratio = 2.7
    try:  # FreeBayes
        ao_index = parts[8].split(":").index("AO")
        ro_index = parts[8].split(":").index("RO")
    except ValueError:
        ao_index, ro_index = None, None
    try:  # VarDict
        af_index = parts[8].split(":").index("AF")
    except ValueError:
        af_index = None
    if af_index is None and ao_index is None:
        raise NotImplementedError("Unexpected format annotations: %s" % parts[0])
    def _calc_freq(item):
        try:
            if ao_index is not None and ro_index is not None:
                ao = sum([int(x) for x in item.split(":")[ao_index].split(",")])
                ro = int(item.split(":")[ro_index])
                freq = ao / float(ao + ro)
            elif af_index is not None:
                freq = float(item.split(":")[af_index])
        except (IndexError, ValueError, ZeroDivisionError):
            freq = 0.0
        return freq
    tumor_freq, normal_freq = _calc_freq(parts[9]), _calc_freq(parts[10])
    return normal_freq <= 0.001 or normal_freq <= tumor_freq / thresh_ratio

for line in iter(inVcf):
    line = line.rstrip('\n\r')
    """Call SOMATIC variants from tumor/normal calls, adding REJECT filters and SOMATIC flag.

    Assumes tumor/normal called with tumor first and normal second, as done in bcbio
    implementation.

    Uses MuTect like somatic filter based on implementation in speedseq:
    https://github.com/cc2qe/speedseq/blob/e6729aa2589eca4e3a946f398c1a2bdc15a7300d/bin/speedseq#L62

    Extracts the genotype likelihoods (GLs) from FreeBayes, which are like phred scores
    except not multiplied by 10.0 (https://en.wikipedia.org/wiki/Phred_quality_score).
    For tumors, we retrieve the best likelihood to not be reference (the first GL) and
    for normal, the best likelhood to be reference.

    After calculating the likelihoods, we compare these to thresholds to pass variants
    at tuned sensitivity/precision. Tuning done on DREAM synthetic 3 dataset evaluations.

    We also check that the frequency of the tumor exceeds the frequency of the normal by
    a threshold to avoid calls that are low frequency in both tumor and normal. This supports
    both FreeBayes and VarDict output frequencies.
    """
    # Thresholds are like phred scores, so 3.5 = phred35
    tumor_thresh, normal_thresh = 3.5, 3.5
    if line.startswith("#CHROM"):
        headers = ['##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic event">',
                   ('##FILTER=<ID=REJECT,Description="Not somatic due to normal call frequency '
                    'or phred likelihoods: tumor: %s, normal %s.">')
                    % (int(tumor_thresh * 10), int(normal_thresh * 10))]
        print >> outVcf, "\n".join(headers) + "\n" + line
    elif line.startswith("#"):
        print >> outVcf, line
    else:
        parts = line.split("\t")
        norm_vals = parts[10].split(':')
        if (not all(x in ('.','0') for x in norm_vals)) and _check_lods(parts, tumor_thresh, normal_thresh) and _check_freqs(parts):
            parts[7] = parts[7] + ";SOMATIC"
        else:
            if parts[6] in set([".", "PASS"]):
                parts[6] = "REJECT"
            else:
                parts[6] += ";REJECT"
        line = "\t".join(parts)
        print >> outVcf, line
