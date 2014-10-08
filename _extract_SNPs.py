
import csv, sys
from collections import namedtuple
import xlrd

wb = xlrd.open_workbook('Table S4 2014_specific_snps.xlsx')
ws = wb.sheet_by_name('2014_specific_snps')

SNP = namedtuple('SNP', ('pos', 'orig', 'base'))

snps = []
for row in range(1, ws.nrows):

    pos = int(ws.cell(row, 0).value)
    orig = ws.cell(row, 1).value
    base = ws.cell(row, 2).value

    snps.append(SNP(pos, orig, base))

snps.sort(key=lambda x: x.pos)

matrix = {}
bases = 'ACGT-'

ident = pos = snpi = None
for line in file('1259657_file_s1/alignments/ebov.mafft.fasta'):

    line = line.strip()

    if line[0] == '>':
        ident = line[1:]
        pos = 1
        snpi = 0
        continue

    while snpi < len(snps) and snps[snpi].pos < pos:
        snpi += 1

    while snpi < len(snps) and snps[snpi].pos < pos + len(line):

        snp = snps[snpi]
        base = line[snp.pos - pos]
        if snp.pos not in matrix:
            matrix[snp.pos] = dict([(b, []) for b in bases])
        matrix[snp.pos].setdefault(base, []).append(ident)

        snpi += 1

    pos += len(line)


if True:
    # find groups with 80 vs 81

    groups = set()

    for pos, d in matrix.items():
        for b, l in d.items():
            if len(l) == 80 or len(l) == 81:
                groups.add(tuple(sorted(l))) # for readability

    l = sorted(groups)
    for i in range(len(l)):
        a = set(l[i])
        for j in range(i + 1, len(l)):
            b = set(l[j])
            print('set_%d[%d] vs set_%d[%d] : +(%s) -(%s)' %
                    (i, len(a), j, len(b), str(a-b), str(b-a)))

    sys.exit(0)


w = csv.writer(sys.stdout)
header = ['pos']
for b in bases:
    header += [b + '.2014', b + '.other']

w.writerow(header)

for snp in snps:
    row = [snp.pos]
    for b in bases:
        l = matrix[snp.pos][b]
        n = sum([1 for ident in l if '2014' in ident])
        m = len(l) - n
        s = ''
        if n >0:
            if 'EBOV_2014_KJ660347' in l:
                s += 'K'
            if 'EBOV_2014_G3687' in l:
                s += 'G'
        row += [str(n)+s, m]

    w.writerow(row)

