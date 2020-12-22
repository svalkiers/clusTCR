from .tools import make_pairs_with_distance, createNetwork


class DistancePairs:

    @staticmethod
    def generate(clustering, data):
        """
        Generates the distance pairs of a given clustering and returns a DistancePairs object.
        It does this by calculating the distance between each pair in the clusters

        WARNING: assumes there are no duplicates

        @params
        items_per_cluster: the average size of the clusters that will be generated
        sort_first: if True, the sequences are first sorted on length before generating the clusters and pairs
        """
        cluster_contents = clustering.get_cluster_contents(include_sequences=False)
        pairs = make_pairs_with_distance(data, cluster_contents)
        return DistancePairs(pairs, data)

    def __init__(self, pairs, data):
        self.pairs = pairs
        self.data = data

    def to_dict(self, include_sequences=True):
        """
        Returns a dictionary with all the pairs that have a given edit distance
        For example
            {
                1: [('CAI', 'CAT'), ..],     # CAI and CAT is a pair with distance 1
                2: [('CAI', 'CTA'), ..],     # CAI and CTA is a pair with distance 2
                ...
            }

        include_sequences:
            if True, will include the full sequences as shown in the examples
            if False, will include IDs of the sequences
        """
        d = {}
        for pair in self.pairs:
            distance = pair[0]
            if distance not in d:
                d[distance] = []
            seq1, seq2 = int(pair[1]), int(pair[2])
            if include_sequences:
                seq1, seq2 = self.data[seq1], self.data[seq2]
            d[distance].append((seq1, seq2))
        
        return d
    
    
    def get_weighted_edges(self, distances = None, weight = 2):
        
        if distances is None:
            distances = self.to_dict()
        
        weighted_edges = []
        for dist in distances.keys():
            weighted = dist ** (-weight)
            for pair in distances[dist]:
                edge = (*pair, weighted)
                weighted_edges.append(edge)
        
        return weighted_edges
    

    def percentage_of_distance1_pairs_found(self):
        """
        Calculates the percentage of pairs with distance 1 found
        Accuracy check against a fully accurate and fast hashing method
        """
        pairs_of_distance1 = len([p for p in self.pairs if p[0] == 1])
        total_pairs_of_distance1 = len(createNetwork(self.data))
        return pairs_of_distance1 / total_pairs_of_distance1

