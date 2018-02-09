//
//  reads_compression.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/5/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//


#include "read_compression.h"
#include <iostream>
#define DEBUG false
#define VERIFY false

using namespace std;

char *reference = NULL;
uint8_t snpInRef[MAX_BP_CHR];
uint32_t cumsumP = 0;


/************************
 * Compress the read
 **********************/
uint32_t compress_read(Arithmetic_stream as, read_models models, read_line samLine, uint8_t chr_change){
    int tempF, PosDiff, chrPos, k;
    uint32_t mask;
    uint16_t maskedReadVal;
    // compress read length (assume int)
    for (k=0;k<4;k++) {
        mask = 0xFF<<(k*8);
        maskedReadVal = (uint8_t)(models->read_length & mask)>>(k*8);
        compress_uint8t(as, models->rlength[k], maskedReadVal);
    }
    //printf("read length is: %d\n", models->read_length);
 
    // Compress sam line
    PosDiff = compress_pos(as, models->pos, models->pos_alpha, samLine->pos, chr_change);
    tempF = compress_flag(as, models->flag, samLine->invFlag);
    //tempF = compress_flag(as, models->flag, 0);
    chrPos = compress_edits(as, models, samLine->edits, samLine->cigar, samLine->read, samLine->pos, PosDiff, tempF, &(samLine->cigarFlags));
    
    if (VERIFY) assert(samLine->pos  == chrPos);

    return 1;
}


/***********************
 * Compress the Flag
 **********************/
uint32_t compress_flag(Arithmetic_stream a, stream_model *F, uint16_t flag){
    
    
    // In this case we need to compress the whole flag, althugh the binary information of whether the
    // read is in reverse or not is the most important. Thus, we return the binary info.
    //we use F[0] as there is no context for the flag.
    
    uint16_t x = 0;
    
    x = flag << 11;
    x >>= 15;
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, F[0], flag);
    
    // Update model
    update_model(F[0], flag);
    
    return x;
    
}

/***********************************
 * Compress the Alphabet of Position
 ***********************************/
uint32_t compress_pos_alpha(Arithmetic_stream as, stream_model *PA, uint32_t x){
    
    uint32_t Byte = 0;
    
    // we encode byte per byte i.e. x = [B0 B1 B2 B3]
    
    // Send B0 to the Arithmetic Stream using the alphabet model
    Byte = x >> 24;
    send_value_to_as(as, PA[0], Byte);
    // Update model
    update_model(PA[0], Byte);
    
    // Send B1 to the Arithmetic Stream using the alphabet model
    Byte = (x & 0x00ff0000) >> 16;
    send_value_to_as(as, PA[1], Byte);
    // Update model
    update_model(PA[1], Byte);
    
    // Send B2 to the Arithmetic Stream using the alphabet model
    Byte = (x & 0x0000ff00) >> 8;
    send_value_to_as(as, PA[2], Byte);
    // Update model
    update_model(PA[2], Byte);
    
    // Send B3 to the Arithmetic Stream using the alphabet model
    Byte = (x & 0x000000ff);
    send_value_to_as(as, PA[3], Byte);
    // Update model
    update_model(PA[3], Byte);
    
    return 1;
    
    
}

/***********************
 * Compress the Position
 **********************/
uint32_t compress_pos(Arithmetic_stream as, stream_model *P, stream_model *PA, uint32_t pos, uint8_t chr_change){
    
    static uint32_t prevPos = 0;
    enum {SMALL_STEP = 0, BIG_STEP = 1};
    int32_t x = 0;
    
    // TODO diferent update models for updating -1 and already seen symbols
    // i.e., SMALL_STEP and BIG_STEP
    
    // Check if we are changing chromosomes.
    if (chr_change)
        prevPos = 0;
    
    
    // Compress the position diference (+ 1 to reserve 0 for new symbols)
    x = pos - prevPos + 1;
    
    if (P[0]->alphaExist[x]){
        // Send x to the Arithmetic Stream
        send_value_to_as(as, P[0], P[0]->alphaMap[x]);
        // Update model
        update_model(P[0], P[0]->alphaMap[x]);
    }
    else{
        
        // Send 0 to the Arithmetic Stream
        send_value_to_as(as, P[0], 0);
        
        // Update model
        update_model(P[0], 0);
        
        // Send the new letter to the Arithmetic Stream using the alphabet model
        compress_pos_alpha(as, PA, x);
        
        // Update the statistics of the alphabet for x
        P[0]->alphaExist[x] = 1;
        P[0]->alphaMap[x] = P[0]->alphabetCard; // We reserve the bin 0 for the new symbol flag
        P[0]->alphabet[P[0]->alphabetCard] = x;
        
        // Update model
        update_model(P[0], P[0]->alphabetCard++);
    }
    
    prevPos = pos;
    
    return x;
}

/****************************
 * Compress the match
 *****************************/
uint32_t compress_match(Arithmetic_stream a, stream_model *M, uint8_t match, uint32_t P){
    
    uint32_t ctx = 0;
    static uint8_t  prevM = 0;
    
    
    // Compute Context
    P = (P != 1)? 0:1;
    //prevP = (prevP > READ_LENGTH)? READ_LENGTH:prevP;
    //prevP = (prevP > READ_LENGTH/4)? READ_LENGTH:prevP;
    
    ctx = (P << 1) | prevM;
    
    //ctx = 0;
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, M[ctx], match);
    
    // Update model
    update_model(M[ctx], match);
    
    prevM = match;
    
    return 1;
}

/*************************
 * Compress the snps
 *************************/
uint32_t compress_snps(Arithmetic_stream a, stream_model *S, uint8_t numSnps){
    
    
    // No context is used for the numSnps for the moment.
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, S[0], numSnps);
    
    // Update model
    update_model(S[0], numSnps);
    
    return 1;
    
}


/********************************
 * Compress the indels
 *******************************/
uint32_t compress_indels(Arithmetic_stream a, stream_model *I, uint8_t numIndels){
    
    
    // Nos context is used for the numIndels for the moment.
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, I[0], numIndels);
    
    // Update model
    update_model(I[0], numIndels);
    
    return 1;
    
}

/*******************************
 * Compress the variations
 ********************************/
uint32_t compress_var(Arithmetic_stream a, stream_model *v, uint32_t pos, uint32_t prevPos, uint32_t flag){
    
    uint32_t ctx = 0;
    
    //flag = 0;
    ctx = prevPos << 1 | flag;
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, v[ctx], pos);
    
    // Update model
    update_model(v[ctx], pos);
    
    return 1;
    
}

/*****************************************
 * Compress the chars
 ******************************************/
uint32_t compress_chars(Arithmetic_stream a, stream_model *c, enum BASEPAIR ref, enum BASEPAIR target){
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, c[ref], target);
    
    // Update model
    update_model(c[ref], target);
    
    return 1;
    
}

/*****************************************
 * Compress the edits
 ******************************************/
uint32_t compress_edits(Arithmetic_stream as, read_models rs, char *edits, char *cigar, char *read, uint32_t P, uint32_t deltaP, uint8_t flag, uint8_t* cigarFlags){
    
    int i = 0;
    uint32_t Dels[MAX_READ_LENGTH];
    ins Insers[MAX_READ_LENGTH];
    snp SNPs[MAX_READ_LENGTH];

    char recCigar[MAX_CIGAR_LENGTH];

    char *tmpcigar, *tmpEdits;
    char *origCigar = cigar;

    // pos in the reference
    cumsumP = cumsumP + deltaP - 1;// DeltaP is 1-based
    
    uint32_t prev_pos = 0;
    uint32_t delta = 0;
    
    // Check if read matches reference
    bool matches = true;
    for (uint32_t i = 0; i < rs->read_length; i++) {
      if (read[i] != reference[P-1 + i]) {
        matches = false;
        break;
      }
    }
    if (matches) {
      compress_match(as, rs->match, 1, deltaP);

      reconstructCigar(Dels, Insers, 0, 0, rs->read_length, recCigar);
      *cigarFlags = 0;
      if (strcmp(recCigar, origCigar) == 0) {
        *cigarFlags = 1;
      }
      return cumsumP;
    }

    struct sequence seq;
    init_sequence(&seq, Dels, Insers, SNPs);

    uint32_t dist = edit_sequence(read, &(reference[P-1]), rs->read_length, rs->read_length, seq);

    // The edit distance from the sam file may be wrong.
    if (dist == 0) {
      compress_match(as, rs->match, 1, deltaP);
      return cumsumP;
    }

    uint32_t numIns = seq.n_ins;
    uint32_t numSnps = seq.n_snps;
    uint32_t numDels = seq.n_dels;


    if (DEBUG) {
      printf("snps %d, dels %d, ins %d\n", numSnps, numDels, numIns);
      assert(numSnps + numDels + numIns > 0);
    }

    compress_match(as, rs->match, 0, deltaP);
    
    // Compress the edits
    if ((numDels | numIns) == 0) {
        compress_snps(as, rs->snps, numSnps);
    }
    else{
        compress_snps(as, rs->snps, 0);
        compress_indels(as, rs->indels, numSnps);
        compress_indels(as, rs->indels, numDels);
        compress_indels(as, rs->indels, numIns);
    }

    // Store the positions and Chars in the corresponding vector
    prev_pos = 0;
    
    for (i = 0; i < numDels; i++){
        if (VERIFY) assert(prev_pos < rs->read_length);
        Dels[i] = Dels[i] - prev_pos;
        compress_var(as, rs->var, Dels[i], prev_pos, flag);
        if (DEBUG) printf("Delete at offset %d, prev %d \n", Dels[i], prev_pos);
        prev_pos += Dels[i];
    }

    prev_pos = 0;
    for (i = 0; i < numIns; i++){
        Insers[i].pos = Insers[i].pos - prev_pos;
        compress_var(as, rs->var, Insers[i].pos, prev_pos, flag);
        compress_chars(as, rs->chars, O, Insers[i].targetChar);
        if (DEBUG) printf("Insert %c at offset %d, prev_pos %d\n", basepair2char(Insers[i].targetChar), Insers[i].pos, prev_pos);
        prev_pos += Insers[i].pos;
    }
    

    prev_pos = 0;
    for (i = 0; i < numSnps; i++){
        SNPs[i].pos = SNPs[i].pos - prev_pos;
        delta = compute_delta_to_first_snp(prev_pos + 1, rs->read_length);

        delta = (delta << BITS_DELTA);
        compress_var(as, rs->var, SNPs[i].pos, delta + prev_pos, flag); 
        snpInRef[cumsumP - 1 + prev_pos + SNPs[i].pos] = 1;

        compress_chars(as, rs->chars, SNPs[i].refChar, SNPs[i].targetChar);
        if (DEBUG) printf("Replace %c with %c offset %d, prev_pos = %d\n", basepair2char(SNPs[i].refChar), basepair2char(SNPs[i].targetChar), SNPs[i].pos, prev_pos);
        prev_pos += SNPs[i].pos;
    }

    reconstructCigar(Dels, Insers, numDels, numIns, rs->read_length, recCigar);
    //printf("%s\n", recCigar);
    *cigarFlags = 0;

    //printf("%d\n", strlen(recCigar));
    //printf("%d\n", strlen(origCigar));
    if (strcmp(recCigar, origCigar) == 0) {
        *cigarFlags = 1;
    }
    return cumsumP;
    
}

uint32_t compute_delta_to_first_snp(uint32_t prevPos, uint32_t readLen){
    
    uint32_t deltaOut;
    uint32_t j = 0;
    
    deltaOut = readLen + 2;
    
    for (j=0;j<readLen - prevPos; j++){
        if (snpInRef[cumsumP - 1 + j + prevPos] == 1){
            deltaOut = j;
            break;
        }
    }
    
    return deltaOut;
}

uint32_t compute_num_digits(uint32_t x){

    //Get the number of digits (We assume readLength < 1000)

    if (x < 10)
        return 1;
    else if (x < 100)
        return 2;
    else if (x < 1000)
        return 3;
    else if (x < 10000)
        return 4;
    else if (x < 100000)
        return 5;
    else if (x < 1000000)
        return 6;
    else if (x < 10000000)
        return 7;
    else if (x < 100000000)
        return 8;
    else
        return 9;

}

void absolute_to_relative(uint32_t *Dels, uint32_t numDels, ins *Insers, uint32_t numIns) {
    // convert to relative
    uint32_t prev_pos = 0;
    uint32_t i;
    if (numDels > 0) {
        prev_pos = Dels[i];
    }
    /*
    for (i = 0; i < numDels; i++) {
        printf("Before Dels %d, pos %d\n", i, Dels[i]);
    }
    for (i = 0; i < numIns; i++) {
        printf("Before Ins %d, pos %d\n", i, Insers[i].pos);
    }*/

    for (i = 1; i < numDels; i++) {
        Dels[i] = Dels[i] - prev_pos;
        prev_pos += Dels[i]; 
    }

    if (numIns > 0) {
        prev_pos = Insers[0].pos;
    }

    for (i = 1; i < numIns; i++) {
        Insers[i].pos = Insers[i].pos - prev_pos;
        prev_pos += Insers[i].pos; 
    }

    /*
    for (i = 0; i < numDels; i++) {
        printf("After Dels %d, pos %d\n", i, Dels[i]);
    }
    for (i = 0; i < numIns; i++) {
        printf("After Ins %d, pos %d\n", i, Insers[i].pos);
    }*/
}
