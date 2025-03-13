#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_heatmap.py

Created on Sun Mar 9 18:44:50 2025

This script generates a heatmap from a presence/absence DataFrame, used in the
Entheome Genome Extraction Pipeline (EGEP) to visualize BUSCO ID distribution
across assemblies. It reads a CSV file and saves the heatmap as a PNG.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""

import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def plot_heatmap(heatmap_df_file, output_png, title):
    """Generate a heatmap from a DataFrame and save it as a PNG.

    Creates a heatmap visualizing the presence/absence of BUSCO IDs across assemblies,
    with assemblies as columns and BUSCO IDs as rows. The heatmap is styled with a
    blue color map and saved to a specified PNG file.

    Args:
        heatmap_df_file (str): Path to the CSV file containing the DataFrame, where
                               the index is BUSCO IDs and columns are assembly IDs.
        output_png (str): Path where the heatmap PNG will be saved.
        title (str): Title to display above the heatmap.

    Returns:
        None: Saves the heatmap to the specified output file and prints a success
              message.

    Raises:
        FileNotFoundError: If the DataFrame file cannot be found.
        pd.errors.EmptyDataError: If the DataFrame file is empty or malformed.
        ValueError: If the DataFrame lacks expected structure (e.g., no columns).
    """
    try:
        heatmap_data = pd.read_csv(heatmap_df_file, index_col=0)
        if heatmap_data.empty or heatmap_data.columns.empty:
            raise ValueError("DataFrame is empty or has no columns")
        
        plt.figure(figsize=(40, 20))
        sns.heatmap(heatmap_data, cmap="Blues", cbar=False)
        plt.title(title, fontsize=16)
        plt.xlabel("Assemblies (Species ID)", fontsize=12)
        plt.ylabel("BUSCO IDs", fontsize=12)
        plt.xticks([])  # Remove x-axis tick labels
        plt.yticks([])  # Remove y-axis tick labels
        plt.tight_layout()
        plt.savefig(output_png, format="png", dpi=300)
        plt.close()
        print(f"PASS:\tHeatmap saved as {output_png}")
    except FileNotFoundError:
        print(f"Error: Heatmap DataFrame file not found: {heatmap_df_file}")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: Heatmap DataFrame file is empty or malformed: {heatmap_df_file}")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: Invalid DataFrame structure: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: Unexpected issue while plotting heatmap: {e}")
        sys.exit(1)


if __name__ == "__main__":
    """Command-line entry point for generating heatmaps.

    Expects three arguments: the DataFrame file path, output PNG path, and heatmap title.
    """
    if len(sys.argv) != 4:
        print("Usage: python plot_heatmap.py <heatmap_df_file> <output_png> <title>")
        sys.exit(1)
    heatmap_df_file = sys.argv[1]
    output_png = sys.argv[2]
    title = sys.argv[3]
    plot_heatmap(heatmap_df_file, output_png, title)
