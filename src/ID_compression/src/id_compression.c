//
//  id_compression.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 12/10/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

// Compression of the IDs -- work based on the compression of ID in Samcomp by Mahonney and Bonfiled (2012)

#include <stdio.h>
#include "sam_block.h"

uint8_t decompress_uint8t(Arithmetic_stream as, stream_model model){
    
    // Send the value to the Arithmetic Stream
    uint8_t c = read_value_from_as(as, model);
    
    // Update model
    update_model(model, c);
    
    return c;
    
}


int compress_uint8t(Arithmetic_stream as, stream_model model, uint8_t c){
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(as, model, c);
    
    // Update model
    update_model(model, c);
    
    return 1;
    
}

int compress_rname(Arithmetic_stream as, rname_models models, char *rname){
    
    static char prev_name[1024] = {0};
    static int prevChar = 0;
    
    uint32_t ctr = 0;
    
    if(strcmp(rname, prev_name) == 0){
        
        compress_uint8t(as, models->same_ref[0], 0);
        return 0;
        
    }
    
    else{
        compress_uint8t(as, models->same_ref[0], 1);
        while (*rname) {
            compress_uint8t(as, models->rname[prevChar], *rname);
            prev_name[ctr++] = *rname;
            prevChar = *rname++;
        }
        compress_uint8t(as, models->rname[prevChar], 0);
        prev_name[ctr] = 0;
        return 1;
    }
    
}

int compress_mapq(Arithmetic_stream as, mapq_models models, uint8_t mapq){
    
    compress_uint8t(as, models->mapq[0], mapq);
    
    return 0;
}

int decompress_mapq(Arithmetic_stream as, mapq_models models, uint8_t *mapq){
    
    *mapq = decompress_uint8t(as, models->mapq[0]);
    
    return 0;
}

int compress_rnext(Arithmetic_stream as, rnext_models models, char *rnext){
    
    static char prev_name[1024] = {0};
    static int prevChar = 0;
    
    uint32_t ctr = 0;
    
    switch (rnext[0]) {
        case '=':
            compress_uint8t(as, models->same_ref[0], 0);
            return 0;
        case '*':
            compress_uint8t(as, models->same_ref[0], 1);
            return 1;
        default:
            compress_uint8t(as, models->same_ref[0], 2);
            break;
    }
    
    
    while (*rnext) {
        compress_uint8t(as, models->rnext[prevChar], *rnext);
        prev_name[ctr++] = *rnext;
        prevChar = *rnext++;
    }
    compress_uint8t(as, models->rnext[prevChar], 0);
    prev_name[ctr] = 0;
    
    return 1;
}

int decompress_rnext(Arithmetic_stream as, rnext_models models, char *rnext){
    
    static char prev_name[1024] = {0};
    static int prevChar = 0;
    
    uint8_t equal_flag = 0, ch;
    
    uint32_t ctr = 0;
    
    equal_flag = decompress_uint8t(as, models->same_ref[0]);
    
    switch (equal_flag) {
        case 0:
            *rnext = '=', rnext++;
            *rnext = 0;
            return 0;
        case 1:
            *rnext = '*', rnext++;
            *rnext = 0;
            return 1;
        default:
            break;
    }
    
    while ( (ch = decompress_uint8t(as, models->rnext[prevChar])) ) {
        
        *rnext = ch, rnext++;
        prev_name[ctr++] = ch;
        prevChar = ch;
    }
    *rnext = 0;
    
    return 2;
    
}

int compress_pnext(Arithmetic_stream as, pnext_models models, uint32_t pos, int32_t tlen, uint32_t pnext, uint8_t rname_rnextDiff, char* cigar){
    
    uint32_t p = 0, pn = 0, cig_op_len = 0, readLength = 0;
    
    if (rname_rnextDiff) {
        
        compress_uint8t(as, models->raw_pnext[0], (pnext >> 0) & 0xff);
        compress_uint8t(as, models->raw_pnext[1], (pnext >> 8) & 0xff);
        compress_uint8t(as, models->raw_pnext[2], (pnext >> 16) & 0xff);
        compress_uint8t(as, models->raw_pnext[3], (pnext >> 24) & 0xff);
        
        return 0;
    }
    
    if (tlen == 0) {
        if (pnext == pos) {
            compress_uint8t(as, models->assumption[0], 0);
        }
        
        else{
            compress_uint8t(as, models->assumption[0], 1);
            compress_pnext_raw(as, models, pos, pnext);
        }
        
        return 0;
    }
    else if (tlen > 0){
        
        p = pos;
        pn = pnext;
        
    }
    
    else{
        p = pnext;
        pn = pos;
        tlen = -tlen;
    }
    
    while (*cigar) {
        cig_op_len = atoi(cigar);
        cigar += compute_num_digits(cig_op_len);
        if (*cigar == 'M' || *cigar == 'D' || *cigar == 'N') {
            readLength += cig_op_len;
        }
        ++cigar;
    }
    
    if (pn == tlen + p - readLength) {
        compress_uint8t(as, models->assumption[0], 0);
    }
    else{
        compress_uint8t(as, models->assumption[0], 1);
        compress_pnext_raw(as, models, p, pn);
    }
    
        
    return 1;
}


int compress_pnext_raw(Arithmetic_stream as, pnext_models models, uint32_t pos, uint32_t pnext){
    
    uint32_t delta = 0;
    
    
    
    if (pnext == 0) {
        compress_uint8t(as, models->zero[0], 0);
        return 0;
    }
    else{
        compress_uint8t(as, models->zero[0], 1);
        
        if (pnext > pos){
            
            delta = pnext - pos;
            compress_uint8t(as, models->sign[0], 0);
        }
        else {
            delta = pos - pnext;
            compress_uint8t(as, models->sign[0], 1);
        }
        
        compress_uint8t(as, models->diff_pnext[0], (delta >> 0) & 0xff);
        compress_uint8t(as, models->diff_pnext[1], (delta >> 8) & 0xff);
        compress_uint8t(as, models->diff_pnext[2], (delta >> 16) & 0xff);
        compress_uint8t(as, models->diff_pnext[3], (delta >> 24) & 0xff);
    }
    
    return 1;
}


int decompress_pnext(Arithmetic_stream as, pnext_models models, uint32_t pos, int32_t tlen, uint32_t readLength, uint32_t *pnext, uint8_t rname_rnextDiff, char* cigar){
    
    uint32_t delta = 0, comp_flag = 0;
    
    
    
   
    comp_flag = decompress_uint8t(as, models->zero[0]);
    
    if ( comp_flag == 0) {
        *pnext = 0;
        return 0;
    }
    else{
        comp_flag = decompress_uint8t(as, models->sign[0]);
        
        delta |= (decompress_uint8t(as, models->diff_pnext[0]) & 0xff) << 0;
        delta |= (decompress_uint8t(as, models->diff_pnext[1]) & 0xff) << 8;
        delta |= (decompress_uint8t(as, models->diff_pnext[2]) & 0xff) << 16;
        delta |= (decompress_uint8t(as, models->diff_pnext[3]) & 0xff) << 24;
        
        if (comp_flag == 0){
            
            *pnext = delta + pos;
        }
        else {
            
            *pnext = pos - delta;
        }
    }
    
    return 1;
}

int compress_tlen(Arithmetic_stream as, tlen_models models, int32_t tlen){
    
    int32_t delta = 0;
    
    if (tlen == 0) {
        compress_uint8t(as, models->zero[0], 0);
        return 0;
    }
    else{
        compress_uint8t(as, models->zero[0], 1);
        
        if (tlen > 0){
            delta = tlen;
            compress_uint8t(as, models->sign[0], 0);
        }
        else {
            delta = -tlen;
            compress_uint8t(as, models->sign[0], 1);
            
        }
        compress_uint8t(as, models->tlen[0], (delta >> 0) & 0xff);
        compress_uint8t(as, models->tlen[1], (delta >> 8) & 0xff);
        compress_uint8t(as, models->tlen[2], (delta >> 16) & 0xff);
        compress_uint8t(as, models->tlen[3], (delta >> 24) & 0xff);
    }
    
    return 1;
}
int decompress_tlen(Arithmetic_stream as, tlen_models models, int32_t* tlen){
    
    int32_t delta = 0;
    uint32_t decomp_flag = 0;
    
    decomp_flag = decompress_uint8t(as, models->zero[0]);
    
    if ( decomp_flag == 0) {
        *tlen = 0;
        return 0;
    }
    else{
        
        decomp_flag = decompress_uint8t(as, models->sign[0]);
        
        delta |= (decompress_uint8t(as, models->tlen[0]) & 0xff) << 0;
        delta |= (decompress_uint8t(as, models->tlen[1]) & 0xff) << 8;
        delta |= (decompress_uint8t(as, models->tlen[2]) & 0xff) << 16;
        delta |= (decompress_uint8t(as, models->tlen[3]) & 0xff) << 24;
        
        if (decomp_flag == 0)
            *tlen = delta;
        else
            *tlen = -delta;
    }
    
    return 1;
}



int compress_id(Arithmetic_stream as, id_models models, char *id){
    
    static char prev_ID[1024] = {0};
    static uint32_t prev_tokens_ptr[1024] = {0};
    uint8_t token_len = 0, match_len = 0;
    uint32_t i = 0, k = 0, tmp = 0, token_ctr = 0, digit_value = 0, digit_model = 0, prev_digit = 0;
    int delta = 0;
    
    char *id_ptr = id, *id_ptr_tok = NULL;
    
    while (*id_ptr != 0) {
        match_len += (*id_ptr == prev_ID[prev_tokens_ptr[token_ctr] + token_len]), token_len++;
        id_ptr_tok = id_ptr + 1;
        
        // Check if the token is a alphabetic word
        if (isalpha(*id_ptr)) {
            
            while ( isalpha( *id_ptr_tok) ){
                
                // compare with the same token from previous ID
                match_len += (*id_ptr_tok == prev_ID[prev_tokens_ptr[token_ctr] + token_len]), token_len++, id_ptr_tok++;
                
            }
            if (match_len == token_len && !isalpha(prev_ID[prev_tokens_ptr[token_ctr] + token_len])) {
                // The token is the same as last ID
                // Encode a token_type ID_MATCH
                compress_uint8t(as, models->token_type[token_ctr], ID_MATCH);
                
            }
            else {
                // Encode a token type ID_ALPHA, the length of the string and the string
                compress_uint8t(as, models->token_type[token_ctr], ID_ALPHA);
                compress_uint8t(as, models->alpha_len[token_ctr], token_len);
                for (k = 0; k < token_len; k++) {
                    compress_uint8t(as, models->alpha_value[token_ctr], *(id_ptr+k));
                }
            }
            
        }
        // check if the token is a run of zeros
        else if (*id_ptr == '0') {
            
            while ( *id_ptr_tok == '0' ){
                
                // compare with the same token from previous ID
                match_len += ('0' == prev_ID[prev_tokens_ptr[token_ctr] + token_len]), token_len++, id_ptr_tok++;
                
            }
            if (match_len == token_len && prev_ID[prev_tokens_ptr[token_ctr] + token_len] != '0') {
                // The token is the same as last ID
                // Encode a token_type ID_MATCH
                compress_uint8t(as, models->token_type[token_ctr], ID_MATCH);
                
            }
            else {
                // Encode a token type ID_ZEROS and the length of the zeros
                compress_uint8t(as, models->token_type[token_ctr], ID_ZEROS);
                compress_uint8t(as, models->zero_run[token_ctr], token_len);
            }
            
        }
        // Check if the token is a number smaller than (1<<30)
        else if (isdigit(*id_ptr)) {
            
            digit_value = (*id_ptr - '0');
            if (*prev_ID != 0){
                prev_digit = prev_ID[prev_tokens_ptr[token_ctr] + token_len -1] - '0';
            }
            
            if (*prev_ID != 0){
                tmp = 1;
                while (isdigit(prev_ID[prev_tokens_ptr[token_ctr] + tmp])) {
                    prev_digit = prev_digit * 10 + (prev_ID[prev_tokens_ptr[token_ctr] + tmp] - '0');
                    tmp++;
                }
                
            }
            
            while ( isdigit(*id_ptr_tok) && digit_value < (1<<30) ){
                digit_value = digit_value * 10 + (*id_ptr_tok - '0');
                //if (*prev_ID != 0){
                //    prev_digit = prev_digit * 10 + (prev_ID[prev_tokens_ptr[token_ctr] + token_len] - '0');
                //}
                // compare with the same token from previous ID
                match_len += (*id_ptr_tok == prev_ID[prev_tokens_ptr[token_ctr] + token_len]), token_len++, id_ptr_tok++;
                
            }
            if ( match_len == token_len && !isdigit(prev_ID[prev_tokens_ptr[token_ctr] + token_len]) ) {
                // The token is the same as last ID
                // Encode a token_type ID_MATCH
                compress_uint8t(as, models->token_type[token_ctr], ID_MATCH);
                
            }
            else if ( (delta = (digit_value - prev_digit)) < 256 && delta > 0){
                compress_uint8t(as, models->token_type[token_ctr], ID_DELTA);
                compress_uint8t(as, models->delta[token_ctr], delta);
                
            }
            else {
                // Encode a token type ID_DIGIT and the value (byte-based)
                compress_uint8t(as, models->token_type[token_ctr], ID_DIGIT);
                digit_model = (token_ctr << 2);
                compress_uint8t(as, models->integer[digit_model | 0], (digit_value >> 0) & 0xff);
                compress_uint8t(as, models->integer[digit_model | 1], (digit_value >> 8) & 0xff);
                compress_uint8t(as, models->integer[digit_model | 2], (digit_value >> 16) & 0xff);
                compress_uint8t(as, models->integer[digit_model | 3], (digit_value >> 24) & 0xff);
            }
        }
        else {
            
            // compare with the same token from previous ID
            //match_len += (*id_ptr == prev_ID[prev_tokens_ptr[token_ctr]]);
            
            if (match_len == token_len) {
                // The token is the same as last ID
                // Encode a token_type ID_MATCH
                compress_uint8t(as, models->token_type[token_ctr], ID_MATCH);
                
            }
            else {
                // Encode a token type ID_CHAR and the char
                compress_uint8t(as, models->token_type[token_ctr], ID_CHAR);
                compress_uint8t(as, models->chars[token_ctr], *id_ptr);
            }
            
        }
        
        prev_tokens_ptr[token_ctr] = i;
        i += token_len;
        id_ptr = id_ptr_tok;
        match_len = 0;
        token_len = 0;
        token_ctr++;
        
    }
    strcpy(prev_ID, id);
    compress_uint8t(as, models->token_type[token_ctr], ID_END);
    
    return 1;
}

int decompress_rname(Arithmetic_stream as, rname_models models, char *rname){
    
    static char prev_name[1024] = {0};
    static int prevChar = 0;
    
    uint8_t chr_change = 0, ch;
    
    uint32_t ctr = 0;
    
    chr_change = decompress_uint8t(as, models->same_ref[0]);
    
    if (chr_change) {
        
        while ( (ch = decompress_uint8t(as, models->rname[prevChar])) ) {
            
            if (ch == '\n') {
                return -1;
            }
            prev_name[ctr++] = ch;
            prevChar = ch;
            *rname = ch, rname++;
        }
        *rname = '\0';

    }
    
    return chr_change;
    
}




int decompress_id(Arithmetic_stream as, id_models model, char *id){
    
    static char prev_ID[1024] = {0};
    static uint32_t prev_tokens_ptr[1024] = {0};
    static uint32_t prev_tokens_len[1024] = {0};
    //char id[1024] = {0};
    uint8_t token_len = 0;
    uint32_t i = 0, k = 0, token_ctr = 0, digit_value = 0;
    uint32_t delta = 0;
    
    enum token_type tok;
    
    id[0]='\0';
    while ( (tok = uint8t2token(decompress_uint8t(as, model->token_type[token_ctr]))) != ID_END ) {
        
        switch (tok) {
            case ID_MATCH:
                memcpy(id+i, &(prev_ID[prev_tokens_ptr[token_ctr]]), prev_tokens_len[token_ctr]);
                token_len = prev_tokens_len[token_ctr];
                break;
            case ID_ALPHA:
                token_len = decompress_uint8t(as, model->alpha_len[token_ctr]);
                for (k = 0; k < token_len; k++) {
                    id[i+k] = decompress_uint8t(as, model->alpha_value[token_ctr]);
                }
                break;
            case ID_DIGIT:
                digit_value = 0;
                digit_value |= (( decompress_uint8t(as, model->integer[(token_ctr << 2) | 0]) & 0xff ) << 0);
                digit_value |= (( decompress_uint8t(as, model->integer[(token_ctr << 2) | 1]) & 0xff ) << 8);
                digit_value |= (( decompress_uint8t(as, model->integer[(token_ctr << 2) | 2]) & 0xff ) << 16);
                digit_value |= (( decompress_uint8t(as, model->integer[(token_ctr << 2) | 3]) & 0xff ) << 24);
                sprintf(id+i, "%u", digit_value);
                token_len = compute_num_digits(digit_value);
                break;
            case ID_DELTA:
                digit_value = 0;
                delta = decompress_uint8t(as, model->delta[token_ctr]);
                memcpy(id+i, &(prev_ID[prev_tokens_ptr[token_ctr]]), prev_tokens_len[token_ctr]);
                id[i+prev_tokens_len[token_ctr]] = '\0';
                digit_value = atoi(id+i) + delta;
                sprintf(id+i, "%u", digit_value);
                token_len = compute_num_digits(digit_value);
                break;
            case ID_ZEROS:
                token_len = decompress_uint8t(as, model->zero_run[token_ctr]);
                memset(id+i, '0', token_len);
                break;
            case ID_CHAR:
                id[i] = decompress_uint8t(as, model->chars[token_ctr]);
                token_len = 1;
                break;
            default:
                break;
        }
        
        prev_tokens_ptr[token_ctr] = i;
        prev_tokens_len[token_ctr] = token_len;
        i += token_len;
        id[i]='\0';
        token_len = 0;
        token_ctr++;
    }
    id[i]='\0';
    strcpy(prev_ID, id);
    //for (kk=i;kk<=1024;kk++){
    //    prev_ID[kk]='\0';
    //}
    //id[i++] = '\n';
    //putc('@', fs);
    //fwrite(id, i, sizeof(char), fs);
    
    return 1;
}

int compress_aux(Arithmetic_stream as, aux_models models, char **aux_str, uint8_t aux_cnt, aux_block aux)
{
    //We first compress the no. of aux fields
    compress_uint8t(as, models->qAux[0], aux_cnt);
    
    char *ptr, *ptr_data;
    uint8_t mappedChar1, mappedChar2, mappedType, rawType;
    
    uint32_t desc_length;
    
    int32_t i,j,k; char *auxPtr;

    uint8_t it;
    uint16_t numericTagType;
    uint8_t most_common_token; //en el resto definido como 32... hay que cambiar.
    for(it = 0; it<aux_cnt; it++) {
        ptr = aux_str[it];
        
        most_common_token = get_most_common_token(aux->most_common, aux->most_common_size, ptr);
        if (most_common_token != aux->most_common_size)
        {
            //aux field found in most common list, flag = 1.
            compress_uint8t(as, models->most_common_flag[0], 1);
            
            //send token
            compress_uint8t(as, models->most_common_values[0], most_common_token);
            
            //globalCnt++;
            //if(globalCnt>400000)printf("%d\n",globalCnt);
            continue;
        }
        //continue;
        
        //Aux field not found in most common list, flag = 0.
        compress_uint8t(as, models->most_common_flag[0], 0);
        
        numericTagType = preprocessTagType(ptr, &rawType);
        if (numericTagType & NOTFOUNDLUTFLAG) {
            //TagType not in LUT, send LUT flag = 1, then values.
            compress_uint8t(as, models->tagtypeLUTflag[0], 1);
            mappedChar1 = (numericTagType & MASKTAGCHAR1) >> 9;
            mappedChar2 = (numericTagType & MASKTAGCHAR2) >> 3;
            mappedType = numericTagType & MASKTYPE;
            
            compress_uint8t(as,models->tag[0],mappedChar1);
            compress_uint8t(as,models->tag[1],mappedChar2);
            
            if(mappedType == TYPELUTLENGTH) {
                //Unknown type (not in typeLUT, very unlikely, is it even possible?) send flag 1, then ascii value (rawType).
                compress_uint8t(as, models->typeLUTflag[0], 1);
                compress_uint8t(as,models->typeRAW[0],rawType);
            } else {
                //Type in typeLUT, send flag 0, then mapped type value.
                compress_uint8t(as, models->typeLUTflag[0], 0);
                compress_uint8t(as,models->typeLUT[0],mappedType);
            }
        } else {
            //TagType in tagtypeLUT, send LUTflag = 0, then value.
            compress_uint8t(as, models->tagtypeLUTflag[0], 0);
            compress_uint8t(as, models->tagtypeLUT[0], numericTagType & 0xFF);
        }
        
        //We compress the actual data (1st approach: byte per byte)
        ptr_data = ptr + 5;
        
        //Crear enum o algo cuando esta tabla este lista:
        //TIPOS:
        //0: int,
        //1: resto.

        //cambiar a switch
        if(ptr[3]=='i') {
            globalCnt++;
            i = strlen(ptr);
            auxPtr = ptr+5;
            //for(j=5;j<i;j++) printf("%c",ptr[j]);
            k=atoi(auxPtr);

            uint8_t first_byte = (k >> 24) & UINT8_MAX;
            uint8_t second_byte = (k >> 16) & UINT8_MAX;
            uint8_t third_byte = (k >> 8) & UINT8_MAX;
            uint8_t fourth_byte = (k) & UINT8_MAX;
            compress_uint8t(as, models->integers[0], first_byte);
            compress_uint8t(as, models->integers[1], second_byte);
            compress_uint8t(as, models->integers[2], third_byte);
            compress_uint8t(as, models->integers[3], fourth_byte);
        } else {
            desc_length = strlen(ptr_data);
            if(desc_length>=UINT16_MAX) {
              desc_length=UINT16_MAX;
            }
            
            uint8_t first_byte = (desc_length >> 8) & UINT8_MAX;
            uint8_t second_byte = (desc_length) & UINT8_MAX;

            compress_uint8t(as, models->descBytes[0], first_byte);
            compress_uint8t(as, models->descBytes[1], second_byte);
            
            uint8_t buff_cnt = 0;
            while (buff_cnt < UINT16_MAX && *ptr_data!=0) {
                compress_uint8t(as, models->iidBytes[0], *ptr_data);
                ptr_data++;
                buff_cnt++;
            }
        }
        
    }
    
    return 1;
}

int decompress_aux(Arithmetic_stream as, aux_block aux, char* finalLine)
{
    //
    uint8_t aux_fields;
    
    uint8_t it, most_common_flag, most_common_token;
    uint8_t tagtypeLUTflag, typeLUTflag, tagtypeLUTindex, mappedChar1, mappedChar2;
    char tagChar1, tagChar2, typeChar;
    uint16_t desc_length = 0;
    char auxFieldString[MAX_AUX_LENGTH] = {0};
    
    char* finalLine_ptr = finalLine;
    
    aux_models models = aux->models;
    
    aux_fields = decompress_uint8t(as, models->qAux[0]);
    
    for(it = 0; it<aux_fields; it++) {
        
        most_common_flag = decompress_uint8t(as, models->most_common_flag[0]);
        
        if(most_common_flag == 1) {
            //The aux field is in the most common list.
            most_common_token = decompress_uint8t(as, models->most_common_values[0]);
            strcpy(finalLine_ptr,aux->most_common[most_common_token]);
            finalLine_ptr += strlen(aux->most_common[most_common_token]);
            *finalLine_ptr = '\t';
            finalLine_ptr++;
            continue;
        }
        //tag & type dec.
        tagtypeLUTflag = decompress_uint8t(as, models->tagtypeLUTflag[0]);
        
        if(tagtypeLUTflag==1) {
            mappedChar1 = decompress_uint8t(as,models->tag[0]);
            mappedChar2 = decompress_uint8t(as,models->tag[1]);
            tagChar1 = inverseCharMap(mappedChar1);
            tagChar2 = inverseCharMap(mappedChar2);
            
            typeLUTflag = decompress_uint8t(as, models->typeLUTflag[0]);
            if(typeLUTflag==1) {
                typeChar = decompress_uint8t(as,models->typeRAW[0]);
            } else {
                typeChar = typeLUT[decompress_uint8t(as,models->typeLUT[0])];
            }
        } else {
            tagtypeLUTindex = decompress_uint8t(as, models->tagtypeLUT[0]);
            tagChar1 = tagTypeLUT[tagtypeLUTindex][0];
            tagChar2 = tagTypeLUT[tagtypeLUTindex][1];
            typeChar = tagTypeLUT[tagtypeLUTindex][2];
        }
        
        
        //cambiar a switch
        uint8_t sign; int value;
        
        char buffer[MAX_AUX_LENGTH] = {0};
        uint16_t buff_cnt;
        
        if(typeChar == 'i') {
            value = decompress_uint8t(as, models->integers[0]) << 24; //alerta (ver analogo en comp)
            value |= decompress_uint8t(as, models->integers[1]) << 16;
            value |= decompress_uint8t(as, models->integers[2]) << 8;
            value |= decompress_uint8t(as, models->integers[3]);
            sprintf(buffer,"%d",value);
            desc_length = strlen(buffer);
            buffer[desc_length] = '\t';
        } else {
            //value dec.
            desc_length = decompress_uint8t(as, models->descBytes[0]) << 8;
            desc_length |= decompress_uint8t(as, models->descBytes[1]);
            for (buff_cnt=0;buff_cnt<desc_length;buff_cnt++) {
                buffer[buff_cnt] = decompress_uint8t(as, models->iidBytes[0]);
            }
            buffer[buff_cnt] = '\t';
        }
        
        // recomp.
        auxFieldString[0] = tagChar1;
        auxFieldString[1] = tagChar2;
        auxFieldString[2] = ':';
        auxFieldString[3] = typeChar;
        auxFieldString[4] = ':';
        
        char* str_pnt;
        str_pnt = auxFieldString+5;
        strcpy(str_pnt,buffer);
        
        strcpy(finalLine_ptr,auxFieldString);
        
        finalLine_ptr += 5+desc_length+1; // tag:type:desc\t
        
    }
    
    *(finalLine_ptr-1) = 0; //We dont want the last tab
    
    return 1;
    
}




int compress_aux_idoia(Arithmetic_stream as, aux_models models, char **aux_str, uint8_t aux_cnt, aux_block aux)
{
    
    char *ptr, *ptr_data;
    uint8_t mappedChar1, mappedChar2, mappedType, rawType;
    
    uint32_t desc_length;
    
    int i,j,k; char *auxPtr;
    
    uint8_t it;
    uint16_t numericTagType, prev_numericTagType, numericTagType_notfound;
    
    uint8_t most_common_token; //en el resto definido como 32... hay que cambiar.
    
    uint8_t use_mostCommonList;
    
    static char tagTypeLUT_update[MAXLUT][4];
    static uint8_t num_auxTypes = TAGTYPELUTLENGTH;
    static uint8_t firstline = 1;
    
    if (firstline == 1){
        for (i=0;i<TAGTYPELUTLENGTH;i++){
            for (j=0;j<4;j++){
                tagTypeLUT_update[i][j] = tagTypeLUT[i][j];
            }
        }
        firstline = 0;
    }
    
    use_mostCommonList = 1;
    
    prev_numericTagType = 0; //first time we use context 0, which is the one for end of line.
    for(it = 0; it<aux_cnt; it++) {
        ptr = aux_str[it];
        
        if (use_mostCommonList == 1){
            
            // MOST COMMON LIST
            most_common_token = get_most_common_token(aux->most_common, aux->most_common_size, ptr);
            if (most_common_token != aux->most_common_size)
            {
                //aux field found in most common list, flag = 1.
                compress_uint8t(as, models->most_common_flag[0], 1);
                
                //send token
                //compress_uint8t(as, models->most_common_values[0], most_common_token);
                compress_uint8t(as, models->most_common_values_wContext[prev_numericTagType], most_common_token);
                
                
                //globalCnt++;
                //if(globalCnt>400000)printf("%d\n",globalCnt);
                
                // STILL NEED TO CHECK THE numericTagType, for next aux fields
                numericTagType = preprocessTagType_update(ptr, &rawType, tagTypeLUT_update, num_auxTypes);
                // numericTagType = preprocessTagType(ptr, &rawType);
                
                if (numericTagType & NOTFOUNDLUTFLAG){
                    prev_numericTagType = 1;
                    // update num_auxTypes and tagTypeLUT_update
                    if (num_auxTypes < MAXLUT){
                        tagTypeLUT_update[num_auxTypes][0] = ptr[0];
                        tagTypeLUT_update[num_auxTypes][1] = ptr[1];
                        tagTypeLUT_update[num_auxTypes][2] = ptr[3];
                        tagTypeLUT_update[num_auxTypes][3] = 0;
                        prev_numericTagType = num_auxTypes + 2;
                        num_auxTypes = num_auxTypes + 1;
                    }
                }else{
                    prev_numericTagType = numericTagType + 2;
                }
                continue;
            }
            compress_uint8t(as, models->most_common_flag[0], 0);
        }
        
        // Check if the TagType is in the LUT
        //numericTagType = preprocessTagType(ptr, &rawType);
        numericTagType = preprocessTagType_update(ptr, &rawType, tagTypeLUT_update, num_auxTypes);
        if (numericTagType & NOTFOUNDLUTFLAG) {
            //TagType not in LUT, encode pos 1
            numericTagType_notfound = numericTagType;
            numericTagType = 1;
            compress_uint8t(as, models->aux_TagType[prev_numericTagType], numericTagType & 0xFF);
            
            //TagType not in LUT, send LUT flag = 1, then values.
            //compress_uint8t(as, models->tagtypeLUTflag[0], 1);
            
            // Since it is not in LUT, need to send values of TAG and TYPE explicitly
            mappedChar1 = (numericTagType_notfound & MASKTAGCHAR1) >> 9;
            mappedChar2 = (numericTagType_notfound & MASKTAGCHAR2) >> 3;
            mappedType = numericTagType_notfound & MASKTYPE;
            
            compress_uint8t(as,models->tag[0],mappedChar1);
            compress_uint8t(as,models->tag[1],mappedChar2);
            
            if(mappedType == TYPELUTLENGTH) {
                //Unknown type (not in typeLUT, very unlikely, is it even possible?) send flag 1, then ascii value (rawType).
                compress_uint8t(as, models->typeLUTflag[0], 1);
                compress_uint8t(as,models->typeRAW[0],rawType);
            } else {
                //Type in typeLUT, send flag 0, then mapped type value.
                compress_uint8t(as, models->typeLUTflag[0], 0);
                compress_uint8t(as,models->typeLUT[0],mappedType);
            }
            
            // update num_auxTypes and tagTypeLUT_update
            if (num_auxTypes < MAXLUT){
                tagTypeLUT_update[num_auxTypes][0] = ptr[0];
                tagTypeLUT_update[num_auxTypes][1] = ptr[1];
                tagTypeLUT_update[num_auxTypes][2] = ptr[3];
                tagTypeLUT_update[num_auxTypes][3] = 0;
                numericTagType = num_auxTypes + 2;
                num_auxTypes = num_auxTypes + 1;
            }
            
        } else {
            // TagType in LUT
            numericTagType = numericTagType + 2; // Remember 0 for end of line, and 1 for not found
            compress_uint8t(as, models->aux_TagType[prev_numericTagType], numericTagType & 0xFF);
            
            //TagType in tagtypeLUT, send LUTflag = 0, then value.
            // compress_uint8t(as, models->tagtypeLUTflag[0], 0);
            //compress_uint8t(as, models->tagtypeLUT[0], numericTagType & 0xFF);
        }
        prev_numericTagType = numericTagType;
        
        //We compress the actual data (1st approach: byte per byte)
        ptr_data = ptr + 5;
        
        //Crear enum o algo cuando esta tabla este lista:
        //TIPOS:
        //0: int,
        //1: resto.
        
        //cambiar a switch
        if(ptr[3]=='i') {
            globalCnt++;
            i = strlen(ptr);
            auxPtr = ptr+5;
            //for(j=5;j<i;j++) printf("%c",ptr[j]);
            k=atoi(auxPtr);
            
            
            k=atoi(auxPtr);
            
            
            // new from geneComp
            uint8_t first_byte = (k >> 24) & UINT8_MAX;
            uint8_t second_byte = (k >> 16) & UINT8_MAX;
            uint8_t third_byte = (k >> 8) & UINT8_MAX;
            uint8_t fourth_byte = (k) & UINT8_MAX;
            compress_uint8t(as, models->integers_wContext[prev_numericTagType], first_byte);
            compress_uint8t(as, models->integers_wContext[(MAXLUT+2)+prev_numericTagType], second_byte);
            compress_uint8t(as, models->integers_wContext[2*(MAXLUT+2)+prev_numericTagType], third_byte);
            compress_uint8t(as, models->integers_wContext[3*(MAXLUT+2)+prev_numericTagType], fourth_byte);
            
            
            /*
             if(k<0) {
             //compress_uint8t(as, models->sign_integers[0], 1);
             compress_uint8t(as, models->sign_integers_wContext[prev_numericTagType], 1);
             k=-k;
             } else {
             //compress_uint8t(as, models->sign_integers[0], 0);
             compress_uint8t(as, models->sign_integers_wContext[prev_numericTagType], 0);
             }
             //compress_uint8t(as, models->integers[0], k); //ALERTA! k > 255?
             compress_uint8t(as, models->integers_wContext[prev_numericTagType], k); //ALERTA! k > 255?
             */
        } else {
            
            // NEw GeneComp
            desc_length = strlen(ptr_data);
            if(desc_length>=UINT16_MAX) {
                desc_length=UINT16_MAX;
            }
            
            uint8_t first_byte = (desc_length >> 8) & UINT8_MAX;
            uint8_t second_byte = (desc_length) & UINT8_MAX;
            
            compress_uint8t(as, models->descBytes_wContext[prev_numericTagType], first_byte);
            compress_uint8t(as, models->descBytes_wContext[(MAXLUT + 2) + prev_numericTagType], second_byte);
            
            uint8_t buff_cnt = 0;
            while (buff_cnt < UINT16_MAX && *ptr_data!=0) {
                compress_uint8t(as, models->iidBytes_wContext[prev_numericTagType], *ptr_data);
                ptr_data++;
                buff_cnt++;
            }
            
            
            /*
             desc_length = strlen(ptr_data);
             if(desc_length>256) desc_length=256;
             
             //compress_uint8t(as, models->descBytes[0], (uint8_t)desc_length);
             compress_uint8t(as, models->descBytes_wContext[prev_numericTagType], (uint8_t)desc_length);
             
             uint8_t buffer[256] = {0};
             uint8_t buff_cnt = 0;
             while (*ptr_data!=0) {
             //compress_uint8t(as, models->iidBytes[0], *ptr_data);
             compress_uint8t(as, models->iidBytes_wContext[prev_numericTagType], *ptr_data);
             ptr_data++;
             }
             */
        }
        
    }
    
    // Now we need to indicate end of AUX fields
    numericTagType = 0;
    if (use_mostCommonList == 1){
        compress_uint8t(as, models->most_common_flag[0], 0);
    }
    compress_uint8t(as, models->aux_TagType[prev_numericTagType], numericTagType & 0xFF);
    
    return 1;
}



int decompress_aux_idoia(Arithmetic_stream as, aux_block aux, char* finalLine)
{
    //
    uint8_t aux_fields;
    
    uint8_t it, most_common_flag, most_common_token;
    uint8_t tagtypeLUTflag, typeLUTflag, tagtypeLUTindex, mappedChar1, mappedChar2, prev_tagtypeLUTindex;
    char tagChar1, tagChar2, typeChar;
    uint8_t desc_length;
    char auxFieldString[MAX_AUX_LENGTH] = {0};
    
    uint8_t rawType;
    
    char* finalLine_ptr = finalLine;
    
    uint16_t tagtypeLUTindex_16;
    
    aux_models models = aux->models;
    
    
    
    //aux_fields = decompress_uint8t(as, models->qAux[0]);
    
    uint8_t moreAux = 1;
    
    uint8_t use_mostCommonList;
    
    int i,j;
    static char tagTypeLUT_update[MAXLUT][4];
    static uint8_t num_auxTypes = TAGTYPELUTLENGTH;
    static uint8_t firstline = 1;
    
    if (firstline == 1){
        for (i=0;i<TAGTYPELUTLENGTH;i++){
            for (j=0;j<4;j++){
                tagTypeLUT_update[i][j] = tagTypeLUT[i][j];
            }
        }
        firstline = 0;
    }
    
    
    
    use_mostCommonList  = 1;
    
    prev_tagtypeLUTindex = 0;
    
    while(moreAux){
        
        if (use_mostCommonList == 1){
            //First check if in most common list
            most_common_flag = decompress_uint8t(as, models->most_common_flag[0]);
            
            if(most_common_flag == 1) {
                //The aux field is in the most common list.
                //most_common_token = decompress_uint8t(as, models->most_common_values[0]);
                most_common_token = decompress_uint8t(as, models->most_common_values_wContext[prev_tagtypeLUTindex]);
                strcpy(finalLine_ptr,aux->most_common[most_common_token]);
                
                //tagtypeLUTindex_16 = preprocessTagType_update(aux->most_common[most_common_token], &rawType, tagTypeLUT_update, num_auxTypes);
                // tagtypeLUTindex_16 = preprocessTagType(finalLine_ptr, &rawType);
                
                finalLine_ptr += strlen(aux->most_common[most_common_token]);
                *finalLine_ptr = '\t';
                finalLine_ptr++;
                
                
                // STILL NEED TO CHECK THE numericTagType, for next aux fields
                //tagtypeLUTindex = preprocessTagType(aux->most_common[most_common_token], &rawType);
                tagtypeLUTindex_16 = preprocessTagType_update(aux->most_common[most_common_token], &rawType, tagTypeLUT_update, num_auxTypes);
                if (tagtypeLUTindex_16 & NOTFOUNDLUTFLAG){
                    prev_tagtypeLUTindex = 1;
                    // update num_auxTypes and tagTypeLUT_update
                    if (num_auxTypes < MAXLUT){
                        tagTypeLUT_update[num_auxTypes][0] = aux->most_common[most_common_token][0];
                        tagTypeLUT_update[num_auxTypes][1] = aux->most_common[most_common_token][1];
                        tagTypeLUT_update[num_auxTypes][2] = aux->most_common[most_common_token][3];
                        tagTypeLUT_update[num_auxTypes][3] = 0;
                        prev_tagtypeLUTindex = num_auxTypes + 2;
                        num_auxTypes = num_auxTypes + 1;
                    }
                }else{
                    prev_tagtypeLUTindex = tagtypeLUTindex_16 + 2;
                }
                
                continue;
            }
        }
        
        
        // Decompress LUTindex + 2
        tagtypeLUTindex = decompress_uint8t(as, models->aux_TagType[prev_tagtypeLUTindex]);
        prev_tagtypeLUTindex = tagtypeLUTindex;
        // Now 3 cases; 0: no more AUX fields 1: TagType not found in LUT else: TagType found in LUT
        if (tagtypeLUTindex == 0){
            moreAux = 0;
            break;
        }else if(tagtypeLUTindex == 1){
            mappedChar1 = decompress_uint8t(as,models->tag[0]);
            mappedChar2 = decompress_uint8t(as,models->tag[1]);
            tagChar1 = inverseCharMap(mappedChar1);
            tagChar2 = inverseCharMap(mappedChar2);
            
            typeLUTflag = decompress_uint8t(as, models->typeLUTflag[0]);
            if(typeLUTflag==1) {
                typeChar = decompress_uint8t(as,models->typeRAW[0]);
            } else {
                typeChar = typeLUT[decompress_uint8t(as,models->typeLUT[0])];
            }
            
            // update num_auxTypes and tagTypeLUT_update
            if (num_auxTypes < MAXLUT){
                tagTypeLUT_update[num_auxTypes][0] = tagChar1;
                tagTypeLUT_update[num_auxTypes][1] = tagChar2;
                tagTypeLUT_update[num_auxTypes][2] = typeChar;
                tagTypeLUT_update[num_auxTypes][3] = 0;
                prev_tagtypeLUTindex = num_auxTypes + 2;
                num_auxTypes = num_auxTypes + 1;
            }
        }else{
            tagtypeLUTindex = tagtypeLUTindex - 2;
            tagChar1 = tagTypeLUT_update[tagtypeLUTindex][0];
            tagChar2 = tagTypeLUT_update[tagtypeLUTindex][1];
            typeChar = tagTypeLUT_update[tagtypeLUTindex][2];
        }
        
        // Decode values if moreAux = 1
        if (moreAux ==1){
            
            //cambiar a switch
            uint8_t sign; int value;
            
            char buffer[MAX_AUX_LENGTH] = {0};
            
            uint16_t buff_cnt;
            
            
            if(typeChar == 'i') {
                
                //NEW
                value = decompress_uint8t(as, models->integers_wContext[prev_tagtypeLUTindex]) << 24; //alerta (ver analogo en comp)
                value |= decompress_uint8t(as, models->integers_wContext[(MAXLUT + 2) + prev_tagtypeLUTindex]) << 16;
                value |= decompress_uint8t(as, models->integers_wContext[2*(MAXLUT + 2) + prev_tagtypeLUTindex]) << 8;
                value |= decompress_uint8t(as, models->integers_wContext[3*(MAXLUT + 2) + prev_tagtypeLUTindex]);
                sprintf(buffer,"%d",value);
                desc_length = strlen(buffer);
                buffer[desc_length] = '\t';
                
                
                /*
                 //sign = decompress_uint8t(as, models->sign_integers[0]);
                 //value = decompress_uint8t(as, models->integers[0]); //alerta (ver analogo en comp)
                 
                 sign = decompress_uint8t(as, models->sign_integers_wContext[prev_tagtypeLUTindex]);
                 value = decompress_uint8t(as, models->integers_wContext[prev_tagtypeLUTindex]); //alerta (ver analogo en comp)
                 
                 if(sign==1) value = -value;
                 sprintf(buffer,"%d",value);
                 desc_length = strlen(buffer);
                 buffer[desc_length] = '\t';
                 */
            } else {
                
                //value dec.
                desc_length = decompress_uint8t(as, models->descBytes_wContext[prev_tagtypeLUTindex]) << 8;
                desc_length |= decompress_uint8t(as, models->descBytes_wContext[(MAXLUT + 2) + prev_tagtypeLUTindex]);
                for (buff_cnt=0;buff_cnt<desc_length;buff_cnt++) {
                    buffer[buff_cnt] = decompress_uint8t(as, models->iidBytes_wContext[prev_tagtypeLUTindex]);
                }
                buffer[buff_cnt] = '\t';
                
                /*
                 //value dec.
                 //desc_length = decompress_uint8t(as, models->descBytes[0]);
                 desc_length = decompress_uint8t(as, models->descBytes_wContext[prev_tagtypeLUTindex]);
                 
                 for (buff_cnt=0;buff_cnt<desc_length;buff_cnt++) {
                 //buffer[buff_cnt] = decompress_uint8t(as, models->iidBytes[0]);
                 buffer[buff_cnt] = decompress_uint8t(as, models->iidBytes_wContext[prev_tagtypeLUTindex]);
                 }
                 buffer[buff_cnt] = '\t';
                 */
            }
            
            // recomp.
            auxFieldString[0] = tagChar1;
            auxFieldString[1] = tagChar2;
            auxFieldString[2] = ':';
            auxFieldString[3] = typeChar;
            auxFieldString[4] = ':';
            
            char* str_pnt;
            str_pnt = auxFieldString+5;
            strcpy(str_pnt,buffer);
            
            strcpy(finalLine_ptr,auxFieldString);
            
            finalLine_ptr += 5+desc_length+1; // tag:type:desc\t
            
        }
        
    }
    
    *(finalLine_ptr-1) = 0; //We dont want the last tab
    
    return 1;
    
}








int compress_cigar(Arithmetic_stream as, read_models models, char *cigar, uint8_t cigarFlags) {
    
    //We check if the cigarFlags == 1. If so, then with the indels is enough to recover the cigar.
    
    //Otherwise, we compress the whole cigar.
    uint8_t cigar_length;
    cigar_length = strlen(cigar);
    
    if(cigarFlags==1) {
        compress_uint8t(as, models->cigarFlags[0], 1);
    } else {
        compress_uint8t(as, models->cigarFlags[0], 0);
        compress_uint8t(as, models->cigar[0],cigar_length);
        while(*cigar) {
            compress_uint8t(as, models->cigar[0],*cigar);
            cigar++;
        }
    }
    
    return 1;
}
