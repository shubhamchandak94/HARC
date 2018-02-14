#ifndef _CLUSTER_H_
#define _CLUSTER_H_

// All of our structures are delcared in lines.h
// Because otherwise it makes this into a definition nightmare

#include "lines.h"

#define MAX_KMEANS_ITERATIONS 1000

// Memory management
struct cluster_list_t *alloc_cluster_list(struct quality_file_t *info);
void free_cluster_list(struct cluster_list_t *);

// Clustering algorithm internals
uint8_t cluster_lines(struct line_block_t *block, struct quality_file_t *info);
double recalculate_means(struct quality_file_t *info);
uint8_t do_cluster_assignment(struct line_t *line, struct quality_file_t *info);
uint8_t assign_cluster(struct line_t *line, struct quality_file_t *info);
void find_distance(struct line_t *line, struct cluster_t *cluster, struct quality_file_t *t);

// Clustering interface
void initialize_kmeans_clustering(struct quality_file_t *info);
void do_kmeans_clustering(struct quality_file_t *info);

#endif
