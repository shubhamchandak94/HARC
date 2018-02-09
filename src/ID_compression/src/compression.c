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
#include "read_compression.h"

//int print_line(struct sam_line_t *sline, uint8_t print_mode, FILE *fs, bool compressing){
//    
//    char foo[] = "CIGAR";
//    
//    int32_t i = 0;
//    
//    switch (print_mode) {
//        case 0:
//            fprintf(fs, "%s", sline->ID);
///*            fprintf(fs, "%d\t", sline->flag);
//            fprintf(fs, "%s\t", sline->rname);
//            fprintf(fs, "%d\t", sline->pos);
//            fprintf(fs, "%d\t", sline->mapq);
//            fprintf(fs, "%s\t", sline->cigar);
//            fprintf(fs, "%s\t", sline->rnext); 
//            fprintf(fs, "%d\t", sline->pnext);
//            fprintf(fs, "%d\t", sline->tlen);
//            fprintf(fs, "%s\t", sline->read);
//            // need to re-reverse quality scores
//            if ((sline->flag & 16) == 16 && !compressing) {
//		//printf("sline print line is: %d\n", sline->readLength);
//                for (i = sline->readLength - 1; i >= 0; --i)
//                    fputc(sline->quals[i], fs);
//                fputc('\t', fs);
//            } else {
//                fprintf(fs, "%s\t", sline->quals);
//            }
//            fprintf(fs, "%s", sline->aux);  
//            fputc('\n', fs); 
//  */          break;
//            
//        default:
//            break;
//    }
//    return 0;
//}
//
//int print_line(struct sam_line_t *sline, uint8_t print_mode, FILE *fs) {
//    return print_line(sline, print_mode, fs, false);
//}

int compress_line(Arithmetic_stream as, sam_block samBlock)  {
    
    // Load the data from the file
    if(load_sam_line(samBlock))
        return 0;
    
    compress_id(as, samBlock->IDs->models, *samBlock->IDs->IDs);
    return 1;
}

//int decompress_line(Arithmetic_stream as, sam_block samBlock, uint8_t lossiness) {
//    
//    int32_t chr_change = 0;
//    
//    uint32_t decompression_flag = 0;
//    
//    struct sam_line_t sline;
//    
//    
//    //This is only for fixed length? i think so.
//  //  sline.readLength = samBlock->read_length;
//    //sline.readLength = samBlock->reads->models->read_length;
//    //printf("sline read length is: %d\n", sline.readLength); 
//    //printf("Decompressing the block...\n");
//    // Loop over the lines of the sam block
//        
////    chr_change = decompress_rname(as, samBlock->rnames->models, sline.rname);
//        
////    if (chr_change == -1)
////        return 0;
//        
//  //  if (chr_change == 1){
//  //          
//        //printf("Chromosome %d decompressed.\n", ++chrCtr);
//            
//        // Store Ref sequence in memory
//  //      store_reference_in_memory(samBlock->fref);
//            
//        // reset cumsumP
// //       cumsumP = 0;
//
//  //      memset(snpInRef, 0, MAX_BP_CHR);
//
//   // }
//
//    decompress_id(as, samBlock->IDs->models, sline.ID);
///*
//    decompress_mapq(as, samBlock->mapq->models, &sline.mapq);
//
//    decompress_rnext(as, samBlock->rnext->models, sline.rnext); 
//
//    decompression_flag = decompress_read(as,samBlock, chr_change, &sline);
//    
//    decompress_cigar(as, samBlock, &sline);
//
//    decompress_tlen(as, samBlock->tlen->models, &sline.tlen);
//
//    decompress_pnext(as, samBlock->pnext->models, sline.pos, sline.tlen, samBlock->read_length, &sline.pnext, sline.rnext[0] != '=', NULL);
//
//    //idoia
//    //decompress_aux(as, samBlock->aux, sline.aux);
//    decompress_aux_idoia(as, samBlock->aux, sline.aux);
//    //idoia
//    
//    if (lossiness == LOSSY) {
//            QVs_decompress(as, samBlock->QVs, decompression_flag, sline.quals);
//    }
//    else
//        QVs_decompress_lossless(as, samBlock->QVs, decompression_flag, sline.quals, (int)strlen(sline.read));
//
//   
//    sline.readLength = samBlock->reads->models->read_length;
//  */  // printf("sline read length before printing is: %d\n", sline.readLength); 
//    print_line(&sline, 0, samBlock->fs);
//
//    return 1;
//}

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


//void* decompress(void *thread_info){
//    
//    uint64_t n = 0;
//    clock_t begin = clock();
//    clock_t ticks;
//    
//
//    struct compressor_info_t *info = (struct compressor_info_t *)thread_info;
//    
//    Arithmetic_stream as = alloc_arithmetic_stream(info->mode, info->fcomp);
//    
//    sam_block samBlock = alloc_sam_models(as, info->fsam, info->fref, info->qv_opts, DECOMPRESSION);
//    
//    //decompress_most_common_list(as, samBlock->aux);
//    
//  //  info->lossiness = decompress_int(as, samBlock->codebook_model);
//    
//    // Start the decompression
//    // initialize the QV model
//  //  if (info->lossiness == LOSSY) {
//   //     initialize_qv_model(as, samBlock->QVs, DECOMPRESSION);
//   // }
//    
//    // Decompress the blocks
//    while(decompress_line(as, samBlock, info->lossiness)){
//        n++;
//    }
//    
//    n += samBlock->block_length;
//    
//    ticks = clock() - begin;
//    printf("Decompression (mapped reads only) took %f\n", ((float)ticks)/CLOCKS_PER_SEC);
//    return NULL;
//}
