from ..modules.prioritcr.analysis import ClusterRepertoire
from ..modules.prioritcr.visualization import ClusterRepertoireVisualization
from ..modules.prioritcr.similarity import tfidf_cosine

from .clustering import ClusteringResult

from functools import partial

def visualize_cluster_graph(
        clusters: ClusteringResult,
        name: str = "cluster_visualization",
        background_color: str = "black"
        ):
    
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
        image_labels = False
        )
    
    # Define distance metric
    # tfidf_cosine requires background repertoire arguments, 
    # which will be passed to the function using partial
    tfidf_cosine_part = partial(
        tfidf_cosine, 
        background_a = cr, 
        background_b = cr
        )
    
    # Create plot
    print("Plotting.")
    if background_color == "black":
        # Black background
        crv.create(
            first_pass_method = "simhash",
            name = name, 
            distance_metric_function = None 
            )
    elif background_color == "white":
        # White background
        crv.create(
            first_pass_method = "simhash",
            name = name, 
            distance_metric_function = None,
            bg_color="#ffffff", 
            line_color="#cccccc" 
            )