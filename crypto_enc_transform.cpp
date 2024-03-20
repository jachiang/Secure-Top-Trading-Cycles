#include "crypto_enc_transform.h"
#include "crypto_utilities.h"

std::vector<Ciphertext<DCRTPoly>> rowToColEnc(std::vector<Ciphertext<DCRTPoly>> &encRows, 
                                              CryptoContext<DCRTPoly> &cryptoContext,
                                              InitRotsMasks &InitRotsMasks,
                                              CryptoOpsLogger &cryptoOpsLogger) {
    // Assumes n x n matrix: n plaintext slots in each row encryption.
    auto n = encRows.size();
    // Generate rotation keys.    
    // std::vector<int32_t> rotIndices;
    // for (size_t i = 0; i < n; i++) { rotIndices.push_back(i); rotIndices.push_back(-i); }
    // cryptoContext->EvalRotateKeyGen(keyPair.secretKey, rotIndices);
    // Populate column containers with encryptions of isolated matrix elements.
    std::vector<std::vector<Ciphertext<DCRTPoly>>> enc_col_container; 
    for (size_t row=0 ; row < n ; ++row){ 
        for (size_t elem=0 ; elem < n ; ++elem){ 
            // Isolate row element and shift element to corresponding position in column.
            // Compute & log Multiplication over ciphertexts.
            auto mask = InitRotsMasks.encMasks();
            TimeVar t; TIC(t);
            auto masked_enc_row = cryptoContext->EvalMult(encRows[row], mask[elem]); // Masked enc(row).
            cryptoContext->ModReduceInPlace(masked_enc_row);
            cryptoOpsLogger.logMult(TOC(t));
            // Compute & log Rotation over ciphertexts.
            TIC(t);
            auto enc_elem = cryptoContext->EvalRotate(masked_enc_row, elem - row);
            cryptoOpsLogger.logRot(TOC(t));
            // Insert isolated column element into column container.
            if (row == 0) {
                std::vector<Ciphertext<DCRTPoly>> enc_elem_vec;
                enc_elem_vec.push_back(enc_elem);
                enc_col_container.push_back(enc_elem_vec); 
            }
            else {
                enc_col_container[elem].push_back(enc_elem);
            }
        }
    }
    // Add all ciphertexts in each column container.
    std::vector<Ciphertext<DCRTPoly>> encCols; 
    for (size_t col=0 ; col < n ; ++col){ 
        // Compute & log AddMany over ciphertexts.
        TimeVar t; TIC(t);
        encCols.push_back(cryptoContext->EvalAddMany(enc_col_container[col]));
        cryptoOpsLogger.logAddMany(TOC(t));
    }   
    return encCols;
}


std::vector<Ciphertext<DCRTPoly>> // Row-encrypted output matrix.
    encElem2Rows(std::vector<std::vector<Ciphertext<DCRTPoly>>> &encMatElems,
                CryptoContext<DCRTPoly> &cryptoContext,
                InitRotsMasks &initRotsMasks,
                CryptoOpsLogger &cryptoOpsLogger) {
    // Derive enc(row) form of input matrix elements.
    auto n = encMatElems.size(); // TODO: Verify input.
    std::vector<Ciphertext<DCRTPoly>> encMatRows;
    for (size_t row=0 ; row < n ; ++row){ 
        std::vector<Ciphertext<DCRTPoly>> encRowContainer;
        for (size_t col=0 ; col < n ; ++col){ 
            // auto encElemMasked = cryptoContext->EvalMult(encMatElems[row][col], initRotsMasks.encMasks()[0]);      
            // cryptoContext->ModReduceInPlace(encElemMasked);
            auto encElemMasked = encMatElems[row][col];
            // Compute & Log Rotation of ciphertexts.
            TimeVar t; TIC(t);
            auto res = cryptoContext->EvalRotate(encElemMasked, -col);
            cryptoOpsLogger.logRot(TOC(t));
            encRowContainer.push_back(res);
        }
        // Compute & Log AddMany over ciphertexts..
        TimeVar t; TIC(t);
        auto res = cryptoContext->EvalAddMany(encRowContainer);
        cryptoOpsLogger.logAddMany(TOC(t));
        encMatRows.push_back(res);
    }
    return encMatRows;       
}
