import sys
import json
import yaml

def isResources(m):
    if 'memory' in m.lower() or 'timeminutes' in m.lower() or 'threads' in m.lower() or 'javaxmx' in m.lower() or 'cores' in m.lower():
        return True
    return False

def isCNVs(m):
    if 'CNVs' in m.lower():
        return True
    return False

def isBWA(m):
    if 'bwa' in m.lower():
        return True
    return False

def isCutadapt(m):
    if 'cutadapt' in m.lower():
        return True
    return False

def isMutect(m):
    if 'mutect2' in m.lower():
        return True
    return False

def isStrelka(m):
    if 'strelka' in m.lower():
        return True
    return False
    
def isVardict(m):
    if 'vardict' in m.lower():
        return True
    return False

def isFastqc(m):
    if 'fastqc' in m.lower():
        return True
    return False

def isMetrics(m):
    if 'metrics' in m.lower():
        return True
    return False

def isBQSR(m):
    if 'bqsr' in m.lower():
        return True
    return False

with open(sys.argv[1]) as inputsf:
    inputs=json.load(inputsf)
    inputkeys = set(inputs)

    with open(sys.argv[2]) as metaf:
        meta=yaml.load(metaf,yaml.Loader)
        metakeys = set(meta['parameter_meta'])
    
    missing=list(inputkeys.difference(metakeys))
    print("missing ({0}):".format(len(missing)))
    missing.sort()
    for m in missing:
        if isResources(m):
            group = 'Resources'
        elif isCNVs(m):
            group = 'CNVs'
        elif isBWA(m):
            group = 'BWAMem'
        elif isCutadapt(m):
            group = 'Cutadapt'
        elif isMutect(m):
            group = 'Mutect'
        elif isStrelka(m):
            group = 'Strelka'
        elif isVardict(m):
            group = 'Vardict'
        elif isFastqc(m):
            group = 'FastQC'
        elif isMetrics(m):
            group = 'Metrics'
        elif isBQSR(m):
            group = 'BQSR'
        else:
            group = 'Other'
        print("\t{}:\n\t\tgroup: {}".format(m,group))
    
    notneeded=list(metakeys.difference(inputkeys))
    print("notneeded ({0}):".format(len(notneeded)))
    notneeded.sort()
    for m in notneeded:
        print("\t{}".format(m))