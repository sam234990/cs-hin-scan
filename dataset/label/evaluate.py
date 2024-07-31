import pandas as pd
import os
import numpy as np
import argparse


def build_community_dict(df_gt):
    community_dict = {}
    for community in df_gt["cluster_id"].unique():
        if community != "o":
            community_dict[community] = set(
                df_gt[df_gt["cluster_id"] == community]["id"].tolist()
            )
    return community_dict


def find_best_f1_score(cluster, df, community_dict):
    best_f1 = 0
    cluster_ids_set = set(df[df["cluster_id"] == cluster]["id"].tolist())

    for community, community_ids_set in community_dict.items():
        TP = len(cluster_ids_set & community_ids_set)

        # Precision & Recall
        precision = TP / len(cluster_ids_set) if len(cluster_ids_set) > 0 else 0
        recall = TP / len(community_ids_set) if len(community_ids_set) > 0 else 0

        # F1 score
        if precision + recall == 0:
            f1 = 0
        else:
            f1 = 2 * (precision * recall) / (precision + recall)

        best_f1 = max(best_f1, f1)

    return best_f1


def calculate_f1_for_o(df, df_gt):
    o_pred_set = set(df[df["cluster_id"] == "o"]["id"].tolist())
    o_gt_set = set(df_gt[df_gt["cluster_id"] == "o"]["id"].tolist())

    TP = len(o_pred_set & o_gt_set)

    # Precision & Recall
    precision = TP / len(o_pred_set) if len(o_pred_set) > 0 else 0
    recall = TP / len(o_gt_set) if len(o_gt_set) > 0 else 0

    # F1 score
    if precision + recall == 0:
        f1 = 0
    else:
        f1 = 2 * (precision * recall) / (precision + recall)

    return f1


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Calculate F1 score for clusters.")
    parser.add_argument(
        "--pred_file", type=str, required=True, help="Path to the prediction file"
    )
    parser.add_argument(
        "--gt_file", type=str, required=True, help="Path to the ground truth file"
    )

    args = parser.parse_args()

    df = pd.read_csv(args.pred_file, sep=" ", header=None, names=["id", "cluster_id"])
    print(df.head(3))
    print(df.shape)

    df_gt = pd.read_csv(args.gt_file, sep=" ", header=None, names=["id", "cluster_id"])
    print(df_gt.shape)

    community_dict = build_community_dict(df_gt)

    clusters = df["cluster_id"].unique()
    print(f"the number of clusters: {len(clusters)}")
    # Remove 'o' from clusters
    clusters = [cluster for cluster in clusters if cluster != "o"]

    f1_scores = {
        cluster: find_best_f1_score(cluster, df, community_dict) for cluster in clusters
    }

    # Calculate F1 score for 'o'
    f1_scores["o"] = calculate_f1_for_o(df, df_gt)

    # Print the top 5 clusters with highest F1 scores
    top_f1_scores = sorted(f1_scores.items(), key=lambda item: item[1], reverse=True)[
        :5
    ]

    print("Top 5 clusters with highest F1 scores:")
    for cluster, score in top_f1_scores:
        print(f"Cluster {cluster}: F1 Score {score}")

    avg_f1 = sum(f1_scores.values()) / len(f1_scores)
    print(avg_f1)
