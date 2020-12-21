def calculate_purity(clustering, epitopes, n):
    contents = clustering.get_cluster_contents()
    purity_amount = 0
    for cluster in contents:
        if len(cluster) <= 1:
            continue
        counts = epitope_counts(cluster, epitopes)
        if len(counts) != 0:
            purity_amount += max(counts.values())
    return purity_amount / n


def calculate_single_cluster_purity(cluster, epitopes):
    return max(epitope_counts(cluster, epitopes).values()) / len(cluster)


def calculate_consistency(clustering, epitopes, n):
    contents = clustering.get_cluster_contents()
    tuples = make_epitope_tuples(contents, epitopes)
    assignments = assign_true_clusters(tuples)
    tcrs_in_true_clusters = sum(assignments.values())
    return tcrs_in_true_clusters / n


def calculate_retention(clustering, epitopes, n):
    contents = clustering.get_cluster_contents()
    count = 0
    for cluster in contents:
        cluster_size = len(cluster)
        if cluster_size > 1:
            count += cluster_size
    return count / n


def make_epitope_tuples(contents, epitopes):
    items = []
    for i, cluster in enumerate(contents):
        items.extend(epitope_tuples(cluster, epitopes, i))
    return sorted(items, key=lambda x: x[1], reverse=True)


def assign_true_clusters(sorted_tuples):
    assignments = dict()
    clusters_assigned = set()

    for item in sorted_tuples:
        epitope, amount, cluster_id = item
        if epitope in assignments or cluster_id in clusters_assigned:
            continue

        assignments[epitope] = amount
        clusters_assigned.add(cluster_id)

    return assignments


def epitope_tuples(cluster, epitopes, cluster_id):
    counts = epitope_counts(cluster, epitopes)
    return [(epitope, amount, cluster_id) for epitope, amount in counts.items()]


def epitope_counts(cluster, epitopes):
    counts = {}
    for id in cluster:
        tcr_epitopes = epitopes[id].split(',')
        for epitope in tcr_epitopes:
            counts[epitope] = counts.get(epitope, 0) + 1
    return counts
