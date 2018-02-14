/**
 * k-means clustering implementation in C
 * 
 * The approach here is parallelizable and should be migrated to a block-based implementation
 * using opencl in order to run faster, or possibly just multithreaded, but I have left it
 * in plain C for the time being to get it working quickly. Note that the means established
 * are discrete values, rather than continuous.
 */

#include "util.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
//#include <malloc.h>

#include "pmf.h"
#include "codebook.h"
#include "cluster.h"

/**
 * Allocate the memory used for the clusters based on the number wanted and column config
 */
struct cluster_list_t *alloc_cluster_list(struct quality_file_t *info) {
	uint8_t j;
	struct cluster_list_t *rtn = (struct cluster_list_t *) calloc(1, sizeof(struct cluster_list_t));

	// Allocate array of cluster structures
	rtn->count = info->cluster_count;
	rtn->clusters = (struct cluster_t *) calloc(info->cluster_count, sizeof(struct cluster_t));
	rtn->distances = (double *) calloc(info->cluster_count, sizeof(double));

	// Fill in each cluster
	for (j = 0; j < info->cluster_count; ++j) {
		rtn->clusters[j].id = j;
		rtn->clusters[j].count = 0;
		rtn->clusters[j].mean = (symbol_t *) calloc(info->columns, sizeof(symbol_t));
		rtn->clusters[j].accumulator = (uint64_t *) calloc(info->columns, sizeof(uint64_t));
		rtn->clusters[j].training_stats = alloc_conditional_pmf_list(info->alphabet, info->columns);
	}

	return rtn;
}

/**
 * Deallocate the memory used for the clusters.
 */
void free_cluster_list(struct cluster_list_t *clusters) {
	uint8_t j;

	for (j = 0; j < clusters->count; ++j) {
		free(clusters->clusters[j].mean);
		free(clusters->clusters[j].accumulator);
		free_conditional_pmf_list(clusters->clusters[j].training_stats);
	}
	free(clusters->distances);
	free(clusters->clusters);
	free(clusters);
}

/**
 * Calculates cluster assignments for the block of lines given, and return
 * status indicating that at least one line changed clusters
 */
uint8_t cluster_lines(struct line_block_t *block, struct quality_file_t *info) {
	uint32_t i;
	uint8_t changed = 0;

	for (i = 0; i < block->count; ++i) {
		changed |= do_cluster_assignment(&block->lines[i], info);
	}

	return changed;
}

/**
 * Updates the cluster means based on their assigned lines. Also clears the line count for
 * the next iteration.
 */
double recalculate_means(struct quality_file_t *info) {
	uint32_t block, line_idx;
	uint32_t i, j;
	struct line_t *line;
	struct cluster_t *cluster;
	uint8_t new_mean;
	double dist, moved;
	double move_max = 0.0;

	// Reset cluster accumulators for new center calculation
	for (i = 0; i < info->cluster_count; ++i) {
		memset(info->clusters->clusters[i].accumulator, 0, info->columns*sizeof(uint64_t));
	}

	// Iterate linewise to accumulate into cluster centers
	for (block = 0; block < info->block_count; ++block) {
		for (line_idx = 0; line_idx < info->blocks[block].count; ++line_idx) {
			line = &info->blocks[block].lines[line_idx];
			cluster = &info->clusters->clusters[line->cluster];
			for (i = 0; i < info->columns; ++i) {
				cluster->accumulator[i] += line->m_data[i];
			}
		}
	}

	// Now find new cluster centers and compute motion
	for (i = 0; i < info->cluster_count; ++i) {
		cluster = &info->clusters->clusters[i];
		dist = 0.0;
		moved = 0.0;

		for (j = 0; j < info->columns; ++j) {
			// Integer division to find the mean, guaranteed to be less than the alphabet size
			new_mean = (uint8_t) (cluster->accumulator[j] / cluster->count);

			// Also figure out how far we've moved
			dist = new_mean - cluster->mean[j];
			moved += dist*dist;

			// Write back the new cluster center
			cluster->mean[j] = new_mean;
		}

		if (moved > move_max)
			move_max = moved;

		if (info->opts->verbose)
			printf("Cluster %d moved %f.\n", i, moved);
	}

	return move_max;
}

/**
 * Compare each line to each cluster to find distances
 */
uint8_t do_cluster_assignment(struct line_t *line, struct quality_file_t *info) {
	uint8_t i;

	for (i = 0; i < info->cluster_count; ++i) {
		find_distance(line, &info->clusters->clusters[i], info);
	}

	return assign_cluster(line, info);
}

/**
 * Assigns a cluster based on the one with the lowest distance
 */
uint8_t assign_cluster(struct line_t *line, struct quality_file_t *info) {
	uint8_t id = 0;
	uint8_t prev_id = line->cluster;
	uint8_t i;
	struct cluster_t *cluster;
	double *distances = info->clusters->distances;
	double d = distances[0];

	// Find the cluster with minimum distance
	for (i = 1; i < info->cluster_count; ++i) {
		if (distances[i] < d) {
			id = i;
			d = distances[i];
		}
	}

	// Assign to that cluster
	line->cluster = id;
	cluster = &info->clusters->clusters[id];
	cluster->count += 1;

	return (prev_id == id) ? 0 : 1;
}

/**
 * Take a line and cluster information and calculates the distance, storing it in the line information vector
 */
void find_distance(struct line_t *line, struct cluster_t *cluster, struct quality_file_t *info) {
	double d = 0.0;
	uint32_t i;
	uint32_t data, mean;

	for (i = 0; i < info->columns; ++i) {
		data = line->m_data[i];
		mean = cluster->mean[i];
		d += (data - mean) * (data - mean);
	}
	info->clusters->distances[cluster->id] = d;
}

/**
 * Initialize the cluster means based on the data given, using random selection
 */
void initialize_kmeans_clustering(struct quality_file_t *info) {
	uint8_t j;
	uint32_t block_id;
	uint32_t line_id;
	struct cluster_list_t *clusters = info->clusters;

	for (j = 0; j < info->cluster_count; ++j) {
		block_id = rand() % info->block_count;
		line_id = rand() % info->blocks[block_id].count;
		memcpy(clusters->clusters[j].mean, info->blocks[block_id].lines[line_id].m_data, info->columns*sizeof(uint8_t));
		if (info->opts->verbose) {
			printf("Chose block %d, line %d.\n", block_id, line_id);
		}
	}
}

/**
 * Do k-means clustering over the set of blocks given to produce a set of clusters that
 * fills the cluster list given
 */
void do_kmeans_clustering(struct quality_file_t *info) {
	uint32_t iter_count = 0;
	uint32_t j;
	uint8_t loop = 1;
	double moved;
	struct cluster_list_t *clusters = info->clusters;

	initialize_kmeans_clustering(info);

	while (iter_count < MAX_KMEANS_ITERATIONS && loop) {
		for (j = 0; j < clusters->count; ++j) {
			clusters->clusters[j].count = 0;
		}

		for (j = 0; j < info->block_count; ++j) {
			cluster_lines(&info->blocks[j], info);
		}

		loop = 0;
		moved = recalculate_means(info);
		if (moved > info->opts->cluster_threshold)
			loop = 1;

		iter_count += 1;
		if (info->opts->verbose) {
			printf("\n");
		}
	}

	if (info->opts->verbose) {
		printf("\nTotal number of iterations: %d.\n", iter_count);
	}
}
