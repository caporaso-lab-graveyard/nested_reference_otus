#!/usr/bin/env python

from cogent.core.tree import PhyloNode

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

def parse_otu_map(maplines):
    """returns {clusterid:[seqids]}"""
    res = []
    for l in maplines:
        fields = l.strip().split('\t')
        res.append((fields[0],fields[1:]))
    return res

def make_nodes(otu_map, length, level):
    """returns nodes that have been created outta the map"""
    nodedist = length / 2.0
    nodes = []
    lookup = {}

    for clusterid, seqids in otu_map:
        rep = seqids[0]
        children = []
        for id_ in seqids:
            node = PhyloNode(Name=id_, Length=nodedist)
            children.append(node)
        parent = PhyloNode(Name="_".join(map(str, [level,clusterid,rep])),
                            Length=None, Children=children)
        nodes.append(parent)
        lookup[rep] = parent
    return lookup,nodes


def join_nodes(parsed):
    """Join the nodes from tips up

    expects parsed to go from high -> low similarity, ie:

    99, 97, 94, ...
    """
    last_lookup, last_nodes = parsed[0]

    for lookup, nodes in parsed[1:]:
        for n in nodes:
            replace_nodes = [last_lookup[c.Name] for c in n.Children]

            new_c = []
            todelete = set([])
            for c in n.Children:
                last_lookup[c.Name].Length = c.Length
                new_c.append(last_lookup[c.Name])
                todelete.add(c)
            
            n.removeDeleted(lambda x: x in todelete)
            
            for c in new_c:
                n.append(c)

        last_lookup = lookup
        last_nodes = nodes

    root = PhyloNode()
    for c in last_nodes:
        root.append(c)

    return root

if __name__ == '__main__':
    from sys import argv

    # expects maps to be in assembly order, ie:
    # gg_99_otu_map.txt,gg_97_otu_map.txt,gg_94_otu_map.txt,...
    otus_in_order = argv[1].split(',')

    last_level = 100
    nodes = []
    for otus in otus_in_order:
        tmp1,level,tmp2,tmp3 = otus.split('_')

        level = int(level)
        length = float(last_level - level)
        otu_map = parse_otu_map(open(otus))
       
        nodes.append(make_nodes(otu_map, length, level))
        
        last_level = level

    tree = join_nodes(nodes)
    f = open(argv[2],'w')
    f.write(tree.getNewick(with_distances=True))
    f.close()

