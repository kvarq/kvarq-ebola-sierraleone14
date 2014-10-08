
VERSION = '0.1'
GENES_COMPATIBILITY = '0.1'

import os.path, functools

from kvarq.genes import Gene, Genotype, Test, Testsuite, Reference, SNP
from kvarq.genes import Genome

# old ebola genome from previous outbreak
EBOV76 = Genome(os.path.join(os.path.dirname(__file__), 'EBOV76.fasta'))


gire14 = Reference('Gire et al (2014) doi 10.1126/science.1259657')

SL1 = Genotype('2014 Sierra Leone outbreak Sublineage 1')
SL2 = Genotype('2014 Sierra Leone outbreak Sublineage 2')
SL3 = Genotype('2014 Sierra Leone outbreak Sublineage 3')

SL1_SNPs = [
        #Test(SNP(genome=EBOV76, pos=127, orig='C', base='T'), SL1, gire14), #80x
        #Test(SNP(genome=EBOV76, pos=155, orig='A', base='C'), SL1, gire14), #80x
        #Test(SNP(genome=EBOV76, pos=182, orig='A', base='G'), SL1, gire14), #80x
        #Test(SNP(genome=EBOV76, pos=187, orig='A', base='G'), SL1, gire14), #80x
        Test(SNP(genome=EBOV76, pos=2263, orig='C', base='T'), SL1, gire14), #80x
        Test(SNP(genome=EBOV76, pos=2314, orig='T', base='C'), SL1, gire14), #80x
        Test(SNP(genome=EBOV76, pos=4340, orig='T', base='C'), SL1, gire14), #80x
        Test(SNP(genome=EBOV76, pos=10057, orig='A', base='G'), SL1, gire14), #80x
        Test(SNP(genome=EBOV76, pos=10065, orig='T', base='G'), SL1, gire14), #80x
        Test(SNP(genome=EBOV76, pos=18764, orig='G', base='A'), SL1, gire14), #80x
    ]

SL2_SNPs = [
        Test(SNP(genome=EBOV76, pos=800, orig='C', base='T'), SL2, gire14), #72x
        Test(SNP(genome=EBOV76, pos=8928, orig='A', base='C'), SL2, gire14), #72x
        Test(SNP(genome=EBOV76, pos=15963, orig='G', base='A'), SL2, gire14), #72x
        Test(SNP(genome=EBOV76, pos=17142, orig='T', base='C'), SL2, gire14), #72x
    ]

SL3_SNPs = [
        Test(SNP(genome=EBOV76, pos=10218, orig='G', base='A'), SL3, gire14), #44x
    ]


class GroupedTestsuite(Testsuite):

    def __init__(self, groups, version):
        tests = functools.reduce(lambda x, y: x + y, groups.values())
        Testsuite.__init__(self, tests, version)
        self.groups = groups

    def __str__(self):
        return 'Showing group(s) of SNPs'

    def _analyse(self, coverages):

        results = []

        for group, tests in self.groups.items():

            n = len(tests)
            x = sum([
                    test.template.validate(coverages[test])
                    for test in tests
                ])
            
            results.append('%s:%d/%d' % (group, x, n))

        return results


sierraleone14 = GroupedTestsuite(dict(SL1=SL1_SNPs, SL2=SL2_SNPs, SL3=SL3_SNPs), VERSION)

