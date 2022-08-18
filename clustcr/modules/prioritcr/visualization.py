import base64
from multiprocessing import Pool
from typing import Callable
import webbrowser
from itertools import repeat

from faerun import Faerun
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Colormap, ListedColormap
import tmap as tm
import faiss
import colorcet as cc


from .analysis import ClusterRepertoire
from .similarity import kmer_jaccard, ClusterRepertoireSimilarity
from .tools import timed


class ClusterRepertoireVisualization:
    def __init__(
        self,
        repertoire: ClusterRepertoire,
        image_labels: bool = False,
        image_labels_method: str = "probability",
    ) -> None:
        self.cr = repertoire
        self.image_labels = image_labels
        self.image_labels_method = image_labels_method
        self._cluster_logo_list = None

    def create(
        self,
        name: str,
        first_pass_method: str,
        first_pass_n: int = 100,
        distance_metric_function: Callable = kmer_jaccard,
        layout_config: tm.LayoutConfiguration = None,
        feature_df: pd.DataFrame = None,
        white_bg: bool = False,
        bg_color: str = None,
        line_color: str = None,
        cmap: Colormap = plt.get_cmap("inferno"),
        simhash_m: int = 128,
    ) -> None:
        """
        Creates and opens an interactive tmap network from the clusters

        Parameters
        ----------
        name : str
            The base filename which will be used to save the network.
        distance_metric_function : Callable, optional
            The cluster distance metric function used, different functions can be found in the
            cluster_similarity module. When using tfidf_cosine, use partial (from functools) to
            add the background repertoire arguments befor passing this function. Default = jaccard
        layout_config : tm.LayoutConfiguration, optional
            additional configuration options. if not passed, defaults will be used. this
            default object can be obtained using .generate_base_layout_config()
        feature_df : pd.DataFrame, optional
            pass this dataframe to customize the features available in-plot.
            a base version can be obtained using the .generate_base_feature_df() method
        image_labels : bool, optional
            whether to include the CDR3 AA logos as images. This can take a while and
            the result is slower to load. default = False
        white_bg : bool, default = False
            plot using a white background (instead of the normal dark grey / black)
        """
        if first_pass_method not in ["minhash", "simhash"]:
            raise ValueError(
                f"{first_pass_method} is not a valid method, choose from exact (pass a distance metric) or simhash"
            )

        if self.image_labels:
            self.x, self.y, self.s, self.t = [None] * 4  # prevents pickling issues
            self._create_cluster_logo_list()

        if layout_config is None:
            layout_config = self.generate_base_layout_config()

        # create layout
        if first_pass_method == "minhash":
            self._construct_network_layout(
                distance_metric_function, first_pass_n, tmap_config=layout_config
            )
        elif first_pass_method == "simhash":
            self._construct_network_layout_simhash(
                distance_metric_function,
                tmap_config=layout_config,
                first_pass_n_nn=first_pass_n,
                m=simhash_m,
            )

        # create plot
        if bg_color is None:
            bg_color = "#f2f2f2" if white_bg else "#111111"

        if line_color is None:
            line_color = "#7F7F7F"

        self._construct_plot(
            name=name,
            feature_df=feature_df,
            image_labels=self.image_labels,
            bg_color=bg_color,
            line_color=line_color,
            cmap=cmap,
        )

    def _construct_network_layout(
        self,
        distance_metric_function: Callable,
        first_pass_n: int,
        tmap_config: tm.LayoutConfiguration,
    ) -> None:

        crs = ClusterRepertoireSimilarity(self.cr)
        edge_list = crs.find_matching_clusters(
            first_pass_n=first_pass_n,
            second_pass_similarity_func=distance_metric_function,
            for_network=True,
        )

        # construct layout
        vertex_count = len(self.cr)
        self.x, self.y, self.s, self.t, _ = tm.layout_from_edge_list(
            vertex_count=vertex_count, edges=edge_list, config=tmap_config
        )

    def _construct_network_layout_simhash(
        self, distance_metric_function, tmap_config, first_pass_n_nn: int = 100, m=128
    ) -> None:
        shash_dict = {c.xid: c._shash(m) for c in self.cr}
        ids = {i: int(c.xid) for i, c in enumerate(self.cr)}

        h = np.array(list(shash_dict.values()))
        arr = np.packbits(np.array(h)).reshape(len(h), m // 8)
        index = faiss.IndexBinaryFlat(m)
        index.add(arr)

        D, I = index.search(arr, first_pass_n_nn)

        edge_list = [
            a
            for b in [
                tuple(zip(repeat(i_1), i_n, score))
                for i_1, i_n, score in zip(range(len(I)), I, D)
            ]
            for a in b
            if a[0] != a[1]
        ]
        if distance_metric_function is None:
            edge_list = [
                (ids[i_1], ids[i_2], score / m) for i_1, i_2, score in edge_list
            ]

        else:
            edge_list = [
                (
                    ids[i_1],
                    ids[i_2],
                    1
                    - distance_metric_function(
                        self.cr.get(ids[i_1]), self.cr.get(ids[i_2])
                    ),
                )
                for i_1, i_2, _ in edge_list
            ]

        # construct plot
        vertex_count = len(self.cr)
        self.x, self.y, self.s, self.t, _ = tm.layout_from_edge_list(
            vertex_count=vertex_count, edges=edge_list, config=tmap_config
        )

    def generate_base_feature_df(self) -> pd.DataFrame:
        df = pd.DataFrame(
            {
                "cluster_id": self._list_id(),
                "normalized_size": self._normalized_size_list(),
                "size": self._list_size(),
                "cdr3_len": [c.sequence_len for c in self.cr],
            }
        ).set_index("cluster_id")
        return df

    def generate_base_layout_config(self) -> tm.LayoutConfiguration:
        return self._default_tmap_config()

    def _list_id(self):
        return [cl.xid for cl in self.cr]

    def _list_size(self):
        return [cl.cluster_size for cl in self.cr]

    def _list_labels(self):
        return [
            f"id: {c.xid}<br>pattern: {c.regex}<br>size: {c.cluster_size}<br>CDR3 length: {c.sequence_len}"
            for c in self.cr
        ]

    def _normalized_size_list(self):
        return np.log10(self._list_size())

    @timed
    def _create_cluster_logo_list(self, n_cpus: int = 8):
        """
        uses multiprocessing to create cluster sequence logos and export them
        to a list of pngs
        ! call this method before constructing network layout, as their resulting
        objects cannot be pickled
        """
        if self._cluster_logo_list is None:
            with Pool(n_cpus) as pool:
                self._cluster_logo_list = pool.map(self._create_cluster_logo, self.cr)
                pool.close()
                pool.join()

    def _create_cluster_logo(self, cluster):
        b = cluster.plot_motif_logo(
            method=self.image_labels_method, export=True, height=0.3, width_per_col=0.15
        )
        img_str = base64.b64encode(b.getvalue())
        return "data:image/bmp;base64," + str(img_str).replace("b'", "").replace(
            "'", ""
        )

    def _default_tmap_config(self) -> tm.LayoutConfiguration:
        cfg = tm.LayoutConfiguration()
        cfg.node_size = 1
        cfg.k = 50
        cfg.kc = 50
        return cfg

    def _is_cat(self, s: pd.Series):
        return pd.api.types.is_categorical_dtype(s) or pd.api.types.is_string_dtype(s)

    def _plot_matplotlib(
        self,
        ax: plt.axes,
        cmap: Colormap,
        line_color: str = "#7F7F7F",
        size_feature: pd.Series = None,
        color_feature: pd.Series = None,
    ):
        coords = {i: v for i, v in enumerate(zip(self.x, self.y))}
        for a, b in list(zip(self.s, self.t)):
            a, b = coords[a], coords[b]
            xc, yc = zip(a, b)
            ax.plot(xc, yc, color=line_color, linewidth=0.5, zorder=0)

        if color_feature is None:
            color_feature = self.generate_base_feature_df()["normalized_size"]
        if size_feature is None:
            size_feature = self.generate_base_feature_df()["normalized_size"] ** 2 * 10

        ax.scatter(
            list(self.x),
            list(self.y),
            s=size_feature,
            c=color_feature,
            cmap=cmap,
        )
        ax.axis("off")
        return ax

    def _parse_metadata(self, feature_df: pd.DataFrame, cmap: Colormap) -> None:
        self.colormap = list()
        self.categorical = list()
        self.legend_title = list()
        self.legend_labels = list()
        self.series_title = list()
        self.c = list()
        self.sizes = list()

        sizes = self._normalized_size_list()
        sizes = sizes / sizes.max() * 6

        for f in feature_df.columns:
            if self._is_cat(feature_df[f]):
                self.categorical.append(True)
                val, lab = feature_df[f].factorize()
                self.c.append(list(val))
                self.legend_labels.append(list(zip(range(len(lab)), lab)))
                self.colormap.append(
                    ListedColormap(
                        sns.color_palette(cc.glasbey_light, as_cmap=True)[: len(lab)]
                    )
                )
            else:
                self.categorical.append(False)
                self.colormap.append(cmap)
                self.c.append(feature_df[f])
                self.legend_labels.append(())
            self.sizes.append(sizes)

            self.legend_title.append(f)
            self.series_title.append(f)

    def _construct_plot(
        self,
        name: str,
        cmap: Colormap,
        feature_df: pd.DataFrame = None,
        image_labels: bool = False,
        bg_color: str = "#111111",
        line_color: str = "#7F7F7F",
    ):

        if image_labels:
            labels = [
                f'id: {c.xid}<br>size: {c.cluster_size}<br>pattern:<br><img src="{im}"/>'
                for c, im in zip(self.cr, self._cluster_logo_list)
            ]
        else:
            labels = self._list_labels()

        if feature_df is None:
            feature_df = self.generate_base_feature_df()
        self._parse_metadata(feature_df=feature_df, cmap=cmap)

        data = dict(x=self.x, y=self.y, labels=labels, s=self.sizes, c=self.c)

        faerun = Faerun(clear_color=bg_color, view="front", coords=False)

        faerun.add_scatter(
            name="tcr_clusters",
            data=data,
            colormap=self.colormap,
            max_point_size=50,
            has_legend=True,
            shader="smoothCircle",
            categorical=self.categorical,
            point_scale=4,
            interactive=True,
            legend_title=self.legend_title,
            legend_labels=self.legend_labels,
            series_title=self.series_title,
        )

        faerun.add_tree(
            name="tcr_clusters_tree",
            data={"from": self.s, "to": self.t},
            point_helper="tcr_clusters",
            color=line_color,
        )

        faerun.plot(file_name=name)
        # webbrowser.open(name + ".html")
