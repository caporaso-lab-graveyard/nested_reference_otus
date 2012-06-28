#1/usr/bin/env python

from cogent.core.tree import PhyloNode
from cogent.util.unit_test import TestCase, main
from unnest import parse_otu_map, make_nodes, join_nodes
from StringIO import StringIO
from cogent.parse.tree import DndParser

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class UnnestTests(TestCase):
    def setUp(self):
        self.clst_99 = parse_otu_map(StringIO(clst_99))
        self.clst_97 = parse_otu_map(StringIO(clst_97))
        self.clst_94 = parse_otu_map(StringIO(clst_94))

    def test_parse_otu_map(self):
        """parse the maps!!"""
        exp = [('0',['10','20','30']),
               ('1',['1','6']),
               ('2',['3']),
               ('3',['8','7'])]
        obs = parse_otu_map(StringIO(clst_99))
        self.assertEqual(obs, exp)

    def test_make_nodes(self):
        """makes nodes..."""
        exp_99_0_10 = PhyloNode(Name="10",Length=0.005)
        exp_99_0_20 = PhyloNode(Name="20",Length=0.005)
        exp_99_0_30 = PhyloNode(Name="30",Length=0.005)
        exp_99_0 = PhyloNode(Name="99_0_10",Length=None, Children=\
                [exp_99_0_10,exp_99_0_20,exp_99_0_30])

        exp_99_1_1 = PhyloNode(Name="1",Length=0.005)
        exp_99_1_6 = PhyloNode(Name="6",Length=0.005)
        exp_99_1 = PhyloNode(Name="99_1_1",Length=None, Children=\
                [exp_99_1_1,exp_99_1_6])

        exp_99_2_3 = PhyloNode(Name="3", Length=0.005)
        exp_99_2 = PhyloNode(Name="99_2_3", Length=None, Children=[exp_99_2_3])

        exp_99_3_8 = PhyloNode(Name="8", Length=0.005)
        exp_99_3_7 = PhyloNode(Name="7", Length=0.005)
        exp_99_3 = PhyloNode(Name="99_3_8", Length=None, Children=\
                [exp_99_3_8,exp_99_3_7])

        exp_lookup = {'10':exp_99_0,'1':exp_99_1,'3':exp_99_2,'8':exp_99_3}

        lookup, nodes = make_nodes(self.clst_99, 0.01, 99)
        self.assertEqual(nodes[0].getNewick(with_distances=True), 
                         exp_99_0.getNewick(with_distances=True))
        self.assertEqual(nodes[1].getNewick(with_distances=True), 
                         exp_99_1.getNewick(with_distances=True))
        self.assertEqual(nodes[2].getNewick(with_distances=True), 
                         exp_99_2.getNewick(with_distances=True))
        self.assertEqual(nodes[3].getNewick(with_distances=True), 
                         exp_99_3.getNewick(with_distances=True))
        self.assertEqual(len(nodes), 4)

        self.assertEqual(lookup.keys(), exp_lookup.keys())
        self.assertEqual(map(str, lookup.values()), 
                         map(str,exp_lookup.values()))

    def test_join_nodes(self):
        """join them nodes! (((99 + 97) + 94) + 91) + ..."""
        parsed = [make_nodes(self.clst_99, 0.01, 99),
                  make_nodes(self.clst_97, 0.02, 97),
                  make_nodes(self.clst_94, 0.03, 94)]

        exp = """((((3:.005)99_2_3:.01,(8:.005,7:.005)99_3_8:.01)97_0_3:.015)94_0_3,
                 (((1:.005,6:.005)99_1_1:.01)97_1_1:.015,
                 ((10:.005,20:.005,30:.005)99_0_10:.01)97_2_10:.015)94_1_1);"""
        expt = DndParser(exp)
        obs = join_nodes(parsed)

        self.assertEqual(obs.getNewick(with_distances=True),
                         expt.getNewick(with_distances=True))
    
clst_99 = """0	10	20	30
1	1	6
2	3
3	8	7
"""
clst_97 = """0	3	8
1	1
2	10
"""
clst_94 = """0	3
1	1	10
"""


if __name__ == '__main__':
    main()
