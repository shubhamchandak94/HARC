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

int print_line(struct sam_line_t *sline, uint8_t print_mode, FILE *fs, bool compressing){
    
    char foo[] = "CIGAR";
    
    int32_t i = 0;
    
    switch (print_mode) {
        case 0:
            fprintf(fs, "%s", sline->ID);
/*            fprintf(fs, "%d\t", sline->flag);
            fprintf(fs, "%s\t", sline->rname);
            fprintf(fs, "%d\t", sline->pos);
            fprintf(fs, "%d\t", sline->mapq);
            fprintf(fs, "%s\t", sline->cigar);
            fprintf(fs, "%s\t", sline->rnext); 
            fprintf(fs, "%d\t", sline->pnext);
            fprintf(fs, "%d\t", sline->tlen);
            fprintf(fs, "%s\t", sline->read);
            // need to re-reverse quality scores
            if ((sline->flag & 16) == 16 && !compressing) {
		//printf("sline print line is: %d\n", sline->readLength);
                for (i = sline->readLength - 1; i >= 0; --i)
                    fputc(sline->quals[i], fs);
                fputc('\t', fs);
            } else {
                fprintf(fs, "%s\t", sline->quals);
            }
            fprintf(fs, "%s", sline->aux);  
            fputc('\n', fs); 
  */          break;
            
        default:
            break;
    }
    return 0;
}

int print_line(struct sam_line_t *sline, uint8_t print_mode, FILE *fs) {
    return print_line(sline, print_mode, fs, false);
}

int compress_line(Arithmetic_stream as, sam_block samBlock, FILE *funmapped, uint8_t lossiness)  {
    
    uint8_t chr_change;
    
    // Load the data from the file
    if(load_sam_line(samBlock))
        return 0;
    //printf("sam line loaded\n");
 
    // If read is unmapped and reference name is *, we assume that all the remaining
    // lines are unmapped and have reference name *.
    // If the read is unmapped but has a position/reference name, we simply use that
    // information to compress it.
/*    if ( (samBlock->reads->lines->invFlag & 4) == 4 && *samBlock->rnames->rnames[0] == '*'){
        read_line line = samBlock->reads->lines; 

        fprintf(funmapped, "%s\t", *samBlock->IDs->IDs);
        fprintf(funmapped, "%d\t", line->invFlag);
        fprintf(funmapped, "%s\t", *samBlock->rnames->rnames);
        fprintf(funmapped, "%d\t", line->pos);
        fprintf(funmapped, "%d\t", *samBlock->mapq->mapq);
        fprintf(funmapped, "%s\t", line->cigar);
        fprintf(funmapped, "%s\t", *samBlock->rnext->rnext); 
        fprintf(funmapped, "%d\t", *samBlock->pnext->pnext);
        fprintf(funmapped, "%d\t", *samBlock->tlen->tlen);
        fprintf(funmapped, "%s\t", line->read);
        int32_t i = 0;
        qv_line_t qline = *samBlock->QVs->qv_lines;
        if ((line->invFlag & 16) == 16) {
            for (i = qline.columns - 1; i >= 0; i--) {
                fputc(qline.data[i] + 33, funmapped);
            }
        } else {
            for (i = 0; i < qline.columns; i++) {
                fputc(qline.data[i] + 33, funmapped);
            }
        }
        fputc('\t', funmapped);
        for (i = 0; i < samBlock->aux->aux_cnt; i++) {
            fprintf(funmapped, "%s", samBlock->aux->aux_str[i]);  
            if (i < samBlock->aux->aux_cnt - 1) fputc('\t', funmapped);
        }
        fputc('\n', funmapped); 
        return 1;
    } else {
        if (unmapped_reads) {
            fprintf(stderr, "compress_line error: There is a mapped read following a read that is not mapped to any chromosome. This probably means that the sam file is not correctly sorted.\n");
            exit(1);
        }
    }
  */      
    // Compress sam line
    
    chr_change = compress_rname(as, samBlock->rnames->models, *samBlock->rnames->rnames);
        
    
    compress_id(as, samBlock->IDs->models, *samBlock->IDs->IDs);
/*
    compress_mapq(as, samBlock->mapq->models, *samBlock->mapq->mapq);

    compress_rnext(as, samBlock->rnext->models, *samBlock->rnext->rnext);

    compress_read(as, samBlock->reads->models, samBlock->reads->lines, chr_change);
    
    compress_cigar(as, samBlock->reads->models, samBlock->reads->lines->cigar, samBlock->reads->lines->cigarFlags);

    compress_tlen(as, samBlock->tlen->models, *samBlock->tlen->tlen);

    compress_pnext_raw(as, samBlock->pnext->models,  samBlock->reads->lines->pos, *samBlock->pnext->pnext);
    
    
    //idoia
    //compress_aux(as, samBlock->aux->models, samBlock->aux->aux_str, samBlock->aux->aux_cnt, samBlock->aux);
    compress_aux_idoia(as, samBlock->aux->models, samBlock->aux->aux_str, samBlock->aux->aux_cnt, samBlock->aux);
    //idoia
    

    if (lossiness == LOSSY)
        QVs_compress(as, samBlock->QVs, samBlock->QVs->qArray);
    else
        QVs_compress_lossless(as, samBlock->QVs->model, samBlock->QVs->qv_lines);
 */   return 1;
}

int decompress_line(Arithmetic_stream as, sam_block samBlock, uint8_t lossiness) {
    
    int32_t chr_change = 0;
    
    uint32_t decompression_flag = 0;
    
    struct sam_line_t sline;
    
    
    //This is only for fixed length? i think so.
  //  sline.readLength = samBlock->read_length;
    //sline.readLength = samBlock->reads->models->read_length;
    //printf("sline read length is: %d\n", sline.readLength); 
    //printf("Decompressing the block...\n");
    // Loop over the lines of the sam block
        
    chr_change = decompress_rname(as, samBlock->rnames->models, sline.rname);
        
    if (chr_change == -1)
        return 0;
        
  //  if (chr_change == 1){
  //          
        //printf("Chromosome %d decompressed.\n", ++chrCtr);
            
        // Store Ref sequence in memory
  //      store_reference_in_memory(samBlock->fref);
            
        // reset cumsumP
 //       cumsumP = 0;

  //      memset(snpInRef, 0, MAX_BP_CHR);

   // }

    decompress_id(as, samBlock->IDs->models, sline.ID);
/*
    decompress_mapq(as, samBlock->mapq->models, &sline.mapq);

    decompress_rnext(as, samBlock->rnext->models, sline.rnext); 

    decompression_flag = decompress_read(as,samBlock, chr_change, &sline);
    
    decompress_cigar(as, samBlock, &sline);

    decompress_tlen(as, samBlock->tlen->models, &sline.tlen);

    decompress_pnext(as, samBlock->pnext->models, sline.pos, sline.tlen, samBlock->read_length, &sline.pnext, sline.rnext[0] != '=', NULL);

    //idoia
    //decompress_aux(as, samBlock->aux, sline.aux);
    decompress_aux_idoia(as, samBlock->aux, sline.aux);
    //idoia
    
    if (lossiness == LOSSY) {
            QVs_decompress(as, samBlock->QVs, decompression_flag, sline.quals);
    }
    else
        QVs_decompress_lossless(as, samBlock->QVs, decompression_flag, sline.quals, (int)strlen(sline.read));

   
    sline.readLength = samBlock->reads->models->read_length;
  */  // printf("sline read length before printing is: %d\n", sline.readLength); 
    print_line(&sline, 0, samBlock->fs);

    return 1;
}

int compress_most_common_list(Arithmetic_stream as, aux_block aux)
{
    uint8_t n,i,l,k;

    aux_models models = aux->models;
    n = aux->most_common_size;
    compress_uint8t(as,models->most_common_list[0],n);

    for (i=0;i<n;i++) {
        l = strlen(aux->most_common[i]);
        compress_uint8t(as,models->most_common_list[0],l);
        for(k=0;k<l;k++)
            compress_uint8t(as,models->most_common_list[0],aux->most_common[i][k]);
    }

    return 1;
}


int decompress_most_common_list(Arithmetic_stream as, aux_block aux)
{
    uint8_t n,i,l,k;

    aux_models models = aux->models;
    n = decompress_uint8t(as,models->most_common_list[0]);

    aux->most_common_size = n;

    char buffer[256];
    for (i=0;i<n;i++) {
        l = decompress_uint8t(as,models->most_common_list[0]);
        for(k=0;k<l;k++)
            buffer[k]=decompress_uint8t(as,models->most_common_list[0]);
        buffer[k]='\0';

        strcpy(aux->most_common[i],buffer);
        //printf("%d -> %s\n",i,aux->most_common[i]);
    }
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
    
    struct qv_options_t opts = *(info.qv_opts);
    
    // Allocs the Arithmetic and the I/O stream
    Arithmetic_stream as = alloc_arithmetic_stream(info.mode, info.fcomp);
    
    // Allocs the different blocks and all the models for the Arithmetic
    sam_block samBlock = alloc_sam_models(as, info.fsam, info.fref, &opts, info.mode);
    
    //idoia
  //  create_most_common_list(samBlock);
   // compress_most_common_list(as, samBlock->aux);
    //idoia
    
 //   if (info.lossiness == LOSSY) {
 //       compress_int(as, samBlock->codebook_model, LOSSY);
 //       initialize_qv_model(as, samBlock->QVs, COMPRESSION);
  //  }
  //  else
   //     compress_int(as, samBlock->codebook_model, LOSSLESS);
   
   printf("start line compression\n"); 
    while (compress_line(as, samBlock, info.funmapped, info.lossiness)) {
        ++lineCtr;
        if (lineCtr % 1000000 == 0) {
          printf("[cbc] compressed %zu lines\n", lineCtr);
        }
    }
    
    // Check if we are in the last block
    compress_rname(as, samBlock->rnames->models, "\n");
    
    //end the compression
    compress_file_size = encoder_last_step(as);
    
    printf("Final Size: %lld\n", compress_file_size);
    
    ticks = clock() - begin;
    
    printf("Compression (mapped reads only) took %f\n", ((float)ticks)/CLOCKS_PER_SEC);
    
    //pthread_exit(NULL);
    return NULL;
}


void* decompress(void *thread_info){
    
    uint64_t n = 0;
    clock_t begin = clock();
    clock_t ticks;
    

    struct compressor_info_t *info = (struct compressor_info_t *)thread_info;
    
    Arithmetic_stream as = alloc_arithmetic_stream(info->mode, info->fcomp);
    
    sam_block samBlock = alloc_sam_models(as, info->fsam, info->fref, info->qv_opts, DECOMPRESSION);
    
    //decompress_most_common_list(as, samBlock->aux);
    
  //  info->lossiness = decompress_int(as, samBlock->codebook_model);
    
    // Start the decompression
    // initialize the QV model
  //  if (info->lossiness == LOSSY) {
   //     initialize_qv_model(as, samBlock->QVs, DECOMPRESSION);
   // }
    
    // Decompress the blocks
    while(decompress_line(as, samBlock, info->lossiness)){
        n++;
    }
    
    n += samBlock->block_length;
    
    ticks = clock() - begin;
    printf("Decompression (mapped reads only) took %f\n", ((float)ticks)/CLOCKS_PER_SEC);
    return NULL;
}
