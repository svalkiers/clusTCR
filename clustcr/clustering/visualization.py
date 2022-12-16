from ..modules.prioritcr.analysis import ClusterRepertoire
from ..modules.prioritcr.visualization import ClusterRepertoireVisualization
from ..modules.prioritcr.similarity import tfidf_cosine, kmer_jaccard, kmer_multiset_jaccard

from .clustering import ClusteringResult

from functools import partial

def visualize_cluster_graph(
        clusters: ClusteringResult,
        name: str = "cluster_visualization",
        background_color: str = "black",
        show_motif_logos: bool = False,
        second_pass_metric: str = None
        ):
    """
    Create interactive visualization of a ClusteringResult.

    Parameters
    ----------
    clusters : ClusteringResult
        ClusTCR ClusteringResult object for which a visualization will be created.
    name : str
        File name of the interactive HTML output. Default is cluster_visualization.
    background_color : str
        Background color of the interactive figure. Can be either black or white. Default is black.
    show_motif_logos : bool
        Construct cluster sequence logos. These will pop up when hovering over the cluster in the interactive HTML.
        Note that construting logos slows down the process. Default is False.
    second_pass_metric : str
        Exhaustive distance calculation of nearby cluster. If set to None, an approximation will be used.
        Available distance metric functions include tfidf_cosine, kmer_jaccard or kmer_multiset_jaccard. 
        Default is None.
    """
    colors = ["black", "white"]
    assert background_color in colors, "Unknown color. Background color must be black or white."
    
    # Calculate cluster similarities
    print("Caclulating cluster similarities.")
    cr = ClusterRepertoire.from_clustcr_result(
        cluster_result = clusters
        )
    # Create interactive visualization background
    print("Initialize interactive visualization.")
    crv = ClusterRepertoireVisualization(
        repertoire = cr,
        image_labels = show_motif_logos
        )

    # Define distance metric, if specified
    if second_pass_metric is not None:
        distance_metric_functions = ["tfidf_cosine", "kmer_jaccard", "kmer_multiset_jaccard"]
        assert second_pass_metric in distance_metric_functions, f"Unknown distance metric '{second_pass_metric}'. Please select one of the following: {distance_metric_functions}."
        if second_pass_metric == "tfidf_cosine":
            # tfidf_cosine requires background repertoire arguments, 
            # which will be passed to the function using partial
            dmf = partial(
                tfidf_cosine, 
                background_a = cr, 
                background_b = cr
                )
        elif second_pass_metric == "kmer_jaccard":
            dmf = kmer_jaccard
        elif second_pass_metric == "kmer_multiset_jaccard":
            dmf = kmer_multiset_jaccard
    else:
        dmf = None
    
    # Create plot
    print("Plotting.")
    if background_color == "black":
        # Black background
        crv.create(
            first_pass_method = "simhash",
            name = name, 
            distance_metric_function = dmf 
            )
    elif background_color == "white":
        # White background
        crv.create(
            first_pass_method = "simhash",
            name = name, 
            distance_metric_function = dmf,
            bg_color="#ffffff", 
            line_color="#cccccc" 
            )