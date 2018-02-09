//
//  compression.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 12/4/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include <stdio.h>
#include <pthread.h>
#include <time.h>
#include <stdbool.h>
#include "sam_block.h"

int print_line(struct sam_line_t *sline, FILE *fs){
    fprintf(fs, "%s", sline->ID);
    return 0;
}

int compress_line(Arithmetic_stream as, sam_block samBlock)  {
    
    // Load the data from the file
    if(load_sam_line(samBlock))
        return 0;
    
    compress_id(as, samBlock->IDs->models, *samBlock->IDs->IDs);
    return 1;
}

int decompress_line(Arithmetic_stream as, sam_block samBlock, FILE *f_id) {
    
    struct sam_line_t sline;
    decompress_id(as, samBlock->IDs->models, sline.ID);
    print_line(&sline, f_id);

    return 1;
}

void* compress(void *thread_info){
    
    uint64_t compress_file_size = 0;
    clock_t begin;
    clock_t ticks;
    
    unsigned long long lineCtr = 0;
    
    printf("Compressing...\n");
    begin = clock();
    
    struct compressor_info_t info = *((struct compressor_info_t *)thread_info);
    
    // Allocs the Arithmetic and the I/O stream
    Arithmetic_stream as = alloc_arithmetic_stream(info.mode, info.fcomp);
    
    // Allocs the different blocks and all the models for the Arithmetic
    sam_block samBlock = alloc_sam_models(as, info.id_array, info.f_order, info.numreads, info.mode);
    
    printf("start line compression\n"); 
    while (compress_line(as, samBlock)) {
        ++lineCtr;
        if (lineCtr % 1000000 == 0) {
          printf("[cbc] compressed %zu lines\n", lineCtr);
        }
    }
    
    //end the compression
    compress_file_size = encoder_last_step(as);
    
    printf("Final Size: %lld\n", compress_file_size);
    
    ticks = clock() - begin;
    
    printf("Compression (mapped reads only) took %f\n", ((float)ticks)/CLOCKS_PER_SEC);
    
    return NULL;
}


void* decompress(void *thread_info){
    
    uint64_t n = 0;
    clock_t begin = clock();
    clock_t ticks;
    

    struct compressor_info_t *info = (struct compressor_info_t *)thread_info;
    
    Arithmetic_stream as = alloc_arithmetic_stream(info->mode, info->fcomp);
    
    sam_block samBlock = alloc_sam_models(as, NULL, NULL, 0, DECOMPRESSION);
    
    // Decompress the blocks
    for(uint32_t n = 0; n < info->numreads; n++) {	
	    decompress_line(as, samBlock, info->f_id);
    }
    
    ticks = clock() - begin;
    printf("Decompression (mapped reads only) took %f\n", ((float)ticks)/CLOCKS_PER_SEC);
    return NULL;
}
