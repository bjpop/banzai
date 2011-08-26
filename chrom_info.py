import sys

from shell_command import shellCommand

def chromInfo(refFile):
    (stdout, stderr, code) = shellCommand('infoseq -noheading ' + refFile)
    if (code != 0):
        print('error from infoseq: ')
        print(stderr)
        sys.exit(code)
    lines = stdout.splitlines()

    # extract the chromosome name and length from the output of infoseq
    chroms = []
    for l in lines:
        words = l.split()
        # just skip bad lines
        if len(words) >= 6:
            chrName = words[2]
            chrLen  = words[5]
            chroms.append((chrName, chrLen))

    return chroms
