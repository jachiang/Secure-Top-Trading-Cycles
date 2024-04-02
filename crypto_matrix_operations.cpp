#include "crypto_matrix_operations.h"

#include <cassert>


std::vector<std::vector<Ciphertext<DCRTPoly>>> // Element-wise-encrypted output matrix
    evalMatrixMul2Pow(std::vector<std::vector<std::vector<Ciphertext<DCRTPoly>>>> &encMatsElems, // Element-wise encrypted input matrices
                      int packingMode, // 0: row/col | 1: full matrix
                      CryptoContext<DCRTPoly> &cryptoContext,
                      InitRotsMasks &initRotsMasks,
                      CryptoOpsLogger &cryptoOpsLogger) {
    
    // Assert number of matrices is power of 2^k for k>0.
    auto numMats = encMatsElems.size();
    int numBits = sizeof(int) * 8; int msbPosition = -1; int bitCtr = 0;
    for (int i = 0; i < numBits; i++) {
        if (numMats >> i & 1) { bitCtr++ ; if (msbPosition == -1){ msbPosition = i + 1; } }}
    // std::cout << bitCtr << " " << msbPosition << " " << numMats << std::endl;
    assert(bitCtr == 1 && msbPosition != 1); // TODO: fails.

    // TODO: Check well-formedness of all input matrices.
    auto n = encMatsElems[0].size();

    // If 2 matrices, multiply.
    if (numMats == 2) {
        // PACKING MODE 0 (ELEMS)
        if (packingMode == 0) {
            std::vector<std::vector<Ciphertext<DCRTPoly>>> encMatElemContainer;
            for (size_t row=0 ; row < n ; ++row){ 
                std::vector<Ciphertext<DCRTPoly>> encMatElemRow;
                for (size_t col=0 ; col < n ; ++col){ 
                    std::vector<Ciphertext<DCRTPoly>> encToAddContainer;
                    for (size_t i=0 ; i < n ; ++i){
                        encToAddContainer.push_back(cryptoContext->EvalMult(encMatsElems[0][row][i], 
                                                                            encMatsElems[1][i][col]));
                    }
                    encMatElemRow.push_back(cryptoContext->EvalAddMany(encToAddContainer));
                }   
                encMatElemContainer.push_back(encMatElemRow);
            }
            return encMatElemContainer;
        }
        // PACKING MODE 1 (ROW/COL)
        else if (packingMode == 1) {
            // TODO: Assert sufficient slots.
            // Convert to row- and col-wise matrix encryptions.
            auto leftEncMat = encElem2Rows(encMatsElems[0],cryptoContext,initRotsMasks,cryptoOpsLogger);
            auto rightEncMat = encElem2Cols(encMatsElems[1],cryptoContext,initRotsMasks,cryptoOpsLogger);

            // Element-wise matrix multiplication.
            std::vector<std::vector<Ciphertext<DCRTPoly>>> encMatElemContainer;
            for (size_t row=0 ; row < n ; ++row){ 
                std::vector<Ciphertext<DCRTPoly>> encMatElemRow;
                for (size_t col=0 ; col < n ; ++col){ 
                    TimeVar t;
                    // Compute & log InnerProduct over ciphertexts.
                    TIC(t); auto encElem = cryptoContext->EvalInnerProduct(leftEncMat[row], rightEncMat[col], n);
                    cryptoOpsLogger.logInnerProd(TOC(t));
                    // Compute & log Multiplication over ciphertexts.
                    TIC(t); auto encElemMasked = cryptoContext->EvalMult(encElem, initRotsMasks.encMasks()[0]);         
                    cryptoContext->ModReduceInPlace(encElemMasked);
                    cryptoOpsLogger.logMult(TOC(t));
                    encMatElemRow.push_back(encElemMasked);
                }
                encMatElemContainer.push_back(encMatElemRow);
            }
            return encMatElemContainer;
        }
        // TODO: PACKING MODE 2 (FULL MATRIX)
        else {
            assert(packingMode == 2);
            // Assert sufficient slots.
            // Dimension of matrix: n
            int k_ceil = std::ceil(std::log2(n));
            int slotsPadded =std::pow(2, k_ceil);

            // Pack left, row-wise encrypted matrix.
            // row1 | row1 | row1 | ... | row2 | row2 | row2 | ... | row3 | row3 | row3 | ... (each row padded 2^{k_ceil})
            std::vector<Ciphertext<DCRTPoly>> encSingleRowContainer;
            for (int row = 0; row < n; row++){
                for (int col = 0; col < n; col++){ 
                    auto pos = row*slotsPadded*slotsPadded+col;
                    encSingleRowContainer.push_back(cryptoContext->EvalRotate(encMatsElems[0][row][col],-pos));
                }
            }
            auto encSingleRows = cryptoContext->EvalAddMany(encSingleRowContainer); 
            auto encRowsReplicated = encSingleRows;
            for (int k = 0; k < k_ceil; k++) {
                auto encRowsReplicated_ = cryptoContext->EvalRotate(encRowsReplicated,-std::pow(2,k)*slotsPadded);
                encRowsReplicated = cryptoContext->EvalAdd(encRowsReplicated,encRowsReplicated_); 
            }
            // Pack right, col-wise encrypted matrix.
            // col1 | col2 | col3 | ... | col1 | col2 | col3 | ... | col1 | col2 | col3 | ... (each row padded 2^{k_ceil})
            std::vector<Ciphertext<DCRTPoly>> encSingleColContainer;
            for (int col = 0; col < n; col++){
                for (int row = 0; row < n; row++){ 
                    auto pos = col*slotsPadded+row;
                    encSingleColContainer.push_back(cryptoContext->EvalRotate(encMatsElems[1][row][col],-pos));
                }
            }
            auto encSingleCols = cryptoContext->EvalAddMany(encSingleColContainer); 
            auto encColsReplicated = encSingleCols;
            for (int k = 0; k < k_ceil; k++) {
                auto encColsReplicated_ = cryptoContext->EvalRotate(encColsReplicated,-std::pow(2,k)*slotsPadded*slotsPadded);
                encColsReplicated = cryptoContext->EvalAdd(encColsReplicated,encColsReplicated_); 
            }
            auto encResMult= cryptoContext->EvalMult(encRowsReplicated, encColsReplicated);
            auto encResInnerProd = evalPrefixAdd(encResMult,slotsPadded,cryptoContext);      

            // Extract individual ciphertexts and rotate to first slot.
            std::vector<std::vector<Ciphertext<DCRTPoly>>> encMatElemContainer;
            for (int row = 0; row < n; row++){
                std::vector<Ciphertext<DCRTPoly>> encElemsRow;
                for (int col = 0; col < n; col++){
                    auto encMaskedElem = cryptoContext->EvalMult(encResInnerProd, initRotsMasks.encMasksFullyPacked()[row*n+col]);
                    cryptoContext->ModReduceInPlace(encMaskedElem);
                    encElemsRow.push_back(cryptoContext->EvalRotate(encMaskedElem,((row*slotsPadded+col)*slotsPadded)));
                }
                encMatElemContainer.push_back(encElemsRow);
            }
            return encMatElemContainer;
        }
    }
    // Otherwise, recursively compute multiplication on left and right half of matrices.
    else {
        std::vector<std::vector<std::vector<Ciphertext<DCRTPoly>>>> encLeftMatsElems(encMatsElems.begin(),encMatsElems.begin()+numMats/2);
        std::vector<std::vector<std::vector<Ciphertext<DCRTPoly>>>> encRightMatsElems(encMatsElems.begin()+numMats/2,encMatsElems.end());
        // std::cout << "Leftmats: " << encLeftMatsElems.size() << std::endl;
        // std::cout << "Rightmats: " << encRightMatsElems.size() << std::endl;
        auto leftMatElems = evalMatrixMul2Pow(encLeftMatsElems, packingMode, cryptoContext, initRotsMasks, cryptoOpsLogger);
        auto rightMatElems = evalMatrixMul2Pow(encRightMatsElems, packingMode, cryptoContext, initRotsMasks, cryptoOpsLogger);
        std::vector<std::vector<std::vector<Ciphertext<DCRTPoly>>>> encMatsContainer;
        encMatsContainer.push_back(leftMatElems); encMatsContainer.push_back(rightMatElems);
        return evalMatrixMul2Pow(encMatsContainer, packingMode, cryptoContext, initRotsMasks, cryptoOpsLogger);
    }
}


std::vector<std::vector<std::vector<Ciphertext<DCRTPoly>>>> // Element-wise-encrypted output matrix.
    evalMatSquarings(std::vector<std::vector<Ciphertext<DCRTPoly>>> &encMatElems,
                    int sqs,
                    int packingMode, // 0: row/col | 1: full matrix
                    CryptoContext<DCRTPoly> &cryptoContext,
                    InitRotsMasks &initRotsMasks,
                    CryptoOpsLogger &cryptoOpsLogger,
                    KeyPair<DCRTPoly> &keyPair // For debugging
                    ) {
    
    // auto n = encMatElems.size();

    std::vector<std::vector<std::vector<Ciphertext<DCRTPoly>>>> encSqMatElems;
    for (int i = 0; i <= sqs; i++) {
        if (i == 0) {
            // Include input matrix as 0'th squaring.
            // std::vector<std::vector<Ciphertext<DCRTPoly>>> encElemContainer;
            // for (size_t row=0 ; row < n ; ++row){ 
            //     std::vector<Ciphertext<DCRTPoly>> encElemRow;
            //     for (size_t col=0 ; col < n ; ++col){ 
            //         // Compute & log Multiplication over ciphertexts.
            //         TimeVar t; 
            //         TIC(t);
            //         auto encElemMasked = cryptoContext->EvalMult(encMatRows[row], initRotsMasks.encMasks()[col]);         
            //         cryptoContext->ModReduceInPlace(encElemMasked);
            //         cryptoOpsLogger.logMult(TOC(t));
            //         // Compute & log Rotation over ciphertexts.
            //         TIC(t); 
            //         auto res = cryptoContext->EvalRotate(encElemMasked, col);
            //         cryptoOpsLogger.logRot(TOC(t));
            //         encElemRow.push_back(res);
            //     }
            //     encElemContainer.push_back(encElemRow);
            // }
            // printEncMatElems(encElemContainer,cryptoContext,keyPair); // TODO: Debugging.
            // encSqMatElems.push_back(encElemContainer);    
            encSqMatElems.push_back(encMatElems);
        }
        // else if (i == 1) {
        //     // First matrix squaring in row-wise encrypted form (due to row-wise enc of input matrix).
        //     // TODO: Use element-wise encrypted input matrix and square.
        //     auto encMatCols = rowToColEnc(encMatRows,cryptoContext,initRotsMasks,cryptoOpsLogger);
        //     std::vector<std::vector<Ciphertext<DCRTPoly>>> encElemContainer;
        //     for (size_t row=0 ; row < n ; ++row){ 
        //         std::vector<Ciphertext<DCRTPoly>> encElemRow;
        //         for (size_t col=0 ; col < n ; ++col){ 
        //             // Compute and log InnerProduct over ciphertexts.
        //             TimeVar t; 
        //             TIC(t);
        //             auto encElem = cryptoContext->EvalInnerProduct(encMatRows[row], encMatCols[col], n);
        //             cryptoOpsLogger.logInnerProd(TOC(t));  
        //             // Compute and log Multiplication over ciphertexts.
        //             TIC(t);
        //             auto encElemMasked = cryptoContext->EvalMult(encElem, initRotsMasks.encMasks()[0]);         
        //             cryptoContext->ModReduceInPlace(encElemMasked);
        //             cryptoOpsLogger.logMult(TOC(t));  
        //             encElemRow.push_back(encElemMasked);
        //         }
        //         encElemContainer.push_back(encElemRow);
        //     }
        //     // printEncMatElems(encElemContainer,cryptoContext,keyPair); // TODO: Debugging.
        //     encSqMatElems.push_back(encElemContainer);
        // }
        else {
            std::vector<std::vector<std::vector<Ciphertext<DCRTPoly>>>> encMatContainer;
            encMatContainer.push_back(encSqMatElems.back());
            encMatContainer.push_back(encSqMatElems.back());
            // TODO: Packing mode. 0:row/col | 1:matrix
            auto encRes = evalMatrixMul2Pow(encMatContainer,packingMode,cryptoContext,initRotsMasks,cryptoOpsLogger);
            // printEncMatElems(encRes,cryptoContext,keyPair); // TODO: Debugging.
            encSqMatElems.push_back(encRes);
        }               
    }
    return encSqMatElems;
}


std::vector<std::vector<Ciphertext<DCRTPoly>>> evalMatSqMul(std::vector<std::vector<Ciphertext<DCRTPoly>>> &encMatElems, 
                                                            int exponent,
                                                            int packingMode, // 0: row/col | 1: full matrix
                                                            CryptoContext<DCRTPoly> &cryptoContext,
                                                            InitRotsMasks &initRotsMasks,
                                                            CryptoOpsLogger &cryptoOpsLogger,
                                                            KeyPair<DCRTPoly> &keyPair // For debugging
                                                            ) {

    // Get msb position of exponent.
    auto numBits = sizeof(int) * 8; int msbPosition = -1; 
    for (int i = 0; i < numBits; i++) { if (exponent >> i & 1) { msbPosition = i+1; } }
    assert(msbPosition > 1); // Required: Exponent > 1.

    // Compute squarings up to msb of exponent. 
    // TODO: Packing mode. 0:row/col | 1:matrix
    auto encMatSqs = evalMatSquarings(encMatElems,msbPosition,packingMode,cryptoContext,initRotsMasks,cryptoOpsLogger,keyPair);

    // Select required squarings.
    std::vector<std::vector<std::vector<Ciphertext<DCRTPoly>>>> encMatSqsActive;
    std::vector<int> activePos;
    for (int i = 0; i < msbPosition; i++) { if (exponent >> i & 1) { encMatSqsActive.push_back(encMatSqs[i]); activePos.push_back(i);} }
    // std::cout << "Selected squarings: " << activePos << std::endl;

    // Multiply selected matrix squarings with log(n) depth.
    // (1) Populate matrix container (encMatsTemp) with multiplications of 2^i matrices.
    std::vector<std::vector<std::vector<Ciphertext<DCRTPoly>>>> encMatsTemp;
    int numSqs = encMatSqsActive.size();
    for (int i = 0; i < numBits; i++) { 
        if (numSqs >> i & 1) {
            std::vector<std::vector<std::vector<Ciphertext<DCRTPoly>>>> encMatsForMult;
            if (i == 0) { encMatsTemp.push_back(encMatSqsActive.back()); encMatSqsActive.pop_back();}
            else {
                for (int j=0; j<std::pow(2,i); j++) {
                    encMatsForMult.push_back(encMatSqsActive.back()); encMatSqsActive.pop_back();
                }
                // TODO: Packing mode. 0:row/col | 1:matrix
                encMatsTemp.push_back(evalMatrixMul2Pow(encMatsForMult,packingMode,cryptoContext,initRotsMasks,cryptoOpsLogger));
            }
        }
    }
    // (2) Sequential multiplication of matrices in container.
    // std::cout << "Number of final mats to multiply sequentially: " << encMatsTemp.size() << std::endl;
    auto encMatRes = encMatsTemp.back(); encMatsTemp.pop_back();
    if (encMatsTemp.size() == 0) { return encMatRes; }
    else {
        for (int i = 0; i < encMatsTemp.size(); i++) {
            // std::cout << "Sequential mult step: " << i << std::endl;
            std::vector<std::vector<std::vector<Ciphertext<DCRTPoly>>>> encMatsForMult;
            encMatsForMult.push_back(encMatRes);
            encMatsForMult.push_back(encMatsTemp.back()); encMatsTemp.pop_back();
            // TODO: Packing mode. 0:row/col | 1:matrix
            encMatRes = evalMatrixMul2Pow(encMatsForMult,packingMode,cryptoContext,initRotsMasks,cryptoOpsLogger);
        }
        return encMatRes;
    }
}


Ciphertext<DCRTPoly> evalMatrixVecMult(std::vector<Ciphertext<DCRTPoly>> &encRows, 
                                       Ciphertext<DCRTPoly> &enc_vec,
                                       CryptoContext<DCRTPoly> &cryptoContext,            
                                       InitRotsMasks &initRotsMasks) {
    // TODO: assert slots in initMatrixVecMult consistent with inputs.
    std::vector<Ciphertext<DCRTPoly>> enc_elements;
    for (size_t row=0 ; row < encRows.size() ; ++row){ 
        auto enc_element = cryptoContext->EvalInnerProduct(encRows[row], enc_vec, 
                                                           encRows.size());
        auto enc_element_masked = cryptoContext->EvalMult(enc_element, initRotsMasks.encMasks()[0]);         
        cryptoContext->ModReduceInPlace(enc_element_masked);
        enc_element = cryptoContext->EvalRotate(enc_element_masked,-row);
        enc_elements.push_back(enc_element);
    }
    return cryptoContext->EvalAddMany(enc_elements);
}


Ciphertext<DCRTPoly> evalVecMatrixMult(Ciphertext<DCRTPoly> &enc_vec,
                                       std::vector<Ciphertext<DCRTPoly>> &encCols, 
                                    //    KeyPair<DCRTPoly> keyPair, // TODO: Only public key required for row to col encryption.
                                       CryptoContext<DCRTPoly> &cryptoContext,            
                                    //    InitMatrixVecMult &initMatrixVecMult,
                                       InitRotsMasks &initRotsMasks,
                                       CryptoOpsLogger &cryptoOpsLogger) {
    // TODO: assert slots in initMatrixVecMult consistent with inputs.
    std::vector<Ciphertext<DCRTPoly>> enc_elements;
    for (size_t col=0 ; col < encCols.size() ; ++col){ 
        // Require col encryptions of matrix.
        // auto encCols = rowToColEnc(encRows,cryptoContext, initRotsMasks,cryptoOpsLogger); 
        auto enc_element = cryptoContext->EvalInnerProduct(encCols[col], enc_vec, 
                                                           encCols.size());
        auto enc_element_masked = cryptoContext->EvalMult(enc_element, initRotsMasks.encMasks()[0]);         
        cryptoContext->ModReduceInPlace(enc_element_masked);
        enc_element = cryptoContext->EvalRotate(enc_element_masked,-col);
        enc_elements.push_back(enc_element);
    }
    return cryptoContext->EvalAddMany(enc_elements);
}


// Matrix exponentiation of n x n matrix, encrypted by rows. Multiplicative depth: log(exponent)
// (Key required to generate rotation keys).
std::vector<Ciphertext<DCRTPoly>> evalMatrixExp(std::vector<Ciphertext<DCRTPoly>> &encRows, 
                                            int exponent,
                                            CryptoContext<DCRTPoly> &cryptoContext,
                                            InitRotsMasks &initRotsMasks,
                                            CryptoOpsLogger &CryptoOpsLogger) {

    auto matrix_dim = encRows.size();
    assert(initRotsMasks.slots == matrix_dim);

    // Generate enc(cols) - rowToColEnc generates rotation keys for context.
    auto encCols = rowToColEnc(encRows, cryptoContext, initRotsMasks, CryptoOpsLogger);

    // Compute all encrypted rows and columns shifted by 0,1,...,n indices.
    std::vector<std::vector<Ciphertext<DCRTPoly>>> encRows_shifted;
    std::vector<std::vector<Ciphertext<DCRTPoly>>> encCols_shifted;
    for (size_t idx=0 ; idx < matrix_dim ; ++idx){ 
        auto enc_row_copy = cryptoContext->EvalRotate(encRows[idx], -matrix_dim); // enc(0...0|row)  
        auto enc_col_copy = cryptoContext->EvalRotate(encCols[idx], -matrix_dim); // enc(0...0|col)
        encRows[idx] = cryptoContext->EvalAdd(encRows[idx], enc_row_copy); // enc(row|row)  
        encCols[idx] = cryptoContext->EvalAdd(encCols[idx], enc_col_copy); // enc(col|col)  
        std::vector<Ciphertext<DCRTPoly>> enc_row_shifted;
        std::vector<Ciphertext<DCRTPoly>> enc_col_shifted;
        for (size_t shift=0 ; shift < matrix_dim ; ++shift){    
            enc_row_shifted.push_back(cryptoContext->EvalRotate(encRows[idx], shift));
            enc_col_shifted.push_back(cryptoContext->EvalRotate(encCols[idx], shift));
        }
        encRows_shifted.push_back(enc_row_shifted);
        encCols_shifted.push_back(enc_col_shifted);
    }
    // n x n Matrix exponentiation. Compute result row by row.
    std::vector<Ciphertext<DCRTPoly>> encRows_res; 
    std::cout << "Matrix exponentation ... " << std::flush; 

        // Iterate through all element indices: for example - i,j,k,l 
        //    (row, i), (row, i+1), (row, i+2) - enc(row),      shifted by i.
        // *  (i, j)  , (i+1, j)  , (i+2, j)   - enc(j'th col), shifted by i. 
        // *  (j, k) ,  (j, k+1)  , (j, k+2)   - enc(j'th row), shifted by k.
        // *  (k, l)  , (k+1, l)  , (k+2, l)   - enc(j'th col), shifted by k. 
        // *  (l, 0) ,  (l, 1)    , (l, 2)     - enc(l'th row).
        std::vector<std::vector<Ciphertext<DCRTPoly>>> enc_add_per_row_container;
        std::vector<Ciphertext<DCRTPoly>> initEncVec;
         for (size_t row=0 ; row < matrix_dim ; ++row){ enc_add_per_row_container.push_back(initEncVec); }
        int indices_dim = exponent-1; VectorIter indices(matrix_dim,indices_dim); bool iterContinue = true;  
        while (iterContinue) { 
            // std::cout << indices.value() << std::endl;
            std::vector<Ciphertext<DCRTPoly>> enc_mult_container; // container for factors of multiplicative terms.
            
            for (int i = 1; i < (indices_dim-2); ++i) { 
                // Push (i, j)  , (i+1, j),   (i+2, j) to container.   
                enc_mult_container.push_back(encCols_shifted[indices.value()[i]][indices.value()[i-1]]); 
                // Push (j, k) ,  (j, k+1)  , (j, k+2) to container.   
                enc_mult_container.push_back(encRows_shifted[indices.value()[i]][indices.value()[i+1]]);
            }
            // Push (k, l)  , (k+1, l)  , (k+2, l) to container.
            enc_mult_container.push_back(encCols_shifted[indices.value().back()][indices.value()[indices_dim-2]]); 
            // Push (l, 0) ,  (l, 1)    , (l, 2) to container.
            enc_mult_container.push_back(encRows_shifted[indices.value().back()][0]);
            // Multiply elements in enc_mult_container, add ciphertext to the enc_add_container.
            auto enc_mult = cryptoContext->EvalMultMany(enc_mult_container);

           // Multiply enc_mult with (row, i), (row, i+1), (row, i+2) ... for all rows.
            for (size_t row=0 ; row < matrix_dim ; ++row){
                enc_add_per_row_container[row].push_back(cryptoContext->EvalMult(enc_mult,
                                                                                 encRows_shifted[row][indices.value()[0]]));
            }

            iterContinue = indices.iterate(); // Continue iterating through i,j,k,l ...
        }

        // Sum all additive terms in enc_add_per_row_container[row], add result to each row encryption.
        for (size_t row=0 ; row < matrix_dim ; ++row){  
            encRows_res.push_back(cryptoContext->EvalAddMany(enc_add_per_row_container[row])); 
        }
        
        std::cout << "completed." << std::endl;
    // }
    return encRows_res;
}

