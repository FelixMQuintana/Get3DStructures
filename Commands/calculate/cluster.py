import numpy as np
import pandas
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage

from Commands.calculate.calculate import CalculateClusters


class HierarchicalClustering(CalculateClusters):
    @staticmethod
    def clustering_algorithm(data: np.ndarray, ordering_file: pandas.DataFrame):
        """

        param data:
        Returns

        """

        dendo = dendrogram(linkage(data), labels=ordering_file[0], color_threshold=6, above_threshold_color='grey',
                           orientation="left")
        plt.axhline(y=6, c='grey', lw=1, linestyle='dashed')