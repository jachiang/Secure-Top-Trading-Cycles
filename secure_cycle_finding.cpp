//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

/*
  Example of a computation circuit of depth 3
  BGVrns demo for a homomorphic multiplication of depth 6 and three different approaches for depth-3 multiplications
 */

#define PROFILE

#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cmath>

#include "openfhe.h"

using namespace lbcrypto;

int modFactorial(int n, int modulus) {
    if (n == 0 || n == 1) { return 1; } 
    else { return ( (n * modFactorial(n - 1, modulus)) % modulus ); }
}

int gcdExtended(int a, int b, int* x, int* y)
{
    if (a == 0) { *x = 0, *y = 1; return b; }
    int x1, y1; int gcd = gcdExtended(b % a, a, &x1, &y1);
    *x = y1 - (b / a) * x1; *y = x1;
    return gcd;
}

int modInverse(int A, int M)
{
    int x, y; int g = gcdExtended(A, M, &x, &y);
    assert(g = 1); // Otherwise, no inverse
    return (x % M + M) % M; 
}

// Helper class for matrix exponentiation: 
// Iterator steps through for all slot indices (i,j,k...) in [0,modulus) x [0,modulus) x [0,modulus) x ...
class vectorIter {
public:
    vectorIter(int modulus, int slots) {
        mod_ = modulus;
        for (size_t i = 0; i < slots; ++i) { vectorIter_.push_back(0); }
    }

    std::vector<int> value() {
        return vectorIter_;
    }

    bool iterate() {
        return iterate_(0);
    }

    bool iterate_(int slot) {
        // Case 0: Highest idx cannot be incremented further.
        if (slot == vectorIter_.size()-1 && (vectorIter_[slot] == mod_-1)) { return false; }
        // Case 1: Increment & carry. 
        else if (vectorIter_[slot] == mod_ - 1) {
            vectorIter_[slot] = 0; if (iterate_(slot+1)) { return true; } else { return false; }}
        // Case 2: Increment.
        else { vectorIter_[slot] = vectorIter_[slot]+1; return true; }
    }

    int mod_;
    std::vector<int> vectorIter_;
};


// Helper method for matrix exponentiation: transforms row encryptions to encryptions of columns.
// Multiplicative depth: 1
// (Key required to generate rotation keys).
std::vector<Ciphertext<DCRTPoly>> rowToColEnc(std::vector<Ciphertext<DCRTPoly>> &enc_rows, 
                                              CryptoContext<DCRTPoly> &cryptoContext,
                                              KeyPair<DCRTPoly> keyPair,
                                              std::vector<Ciphertext<DCRTPoly>> encMasks) {
    // Assumes n x n matrix: n plaintext slots in each row encryption.
    auto n = enc_rows.size();
    // Generate rotation keys.    
    std::vector<int32_t> rotIndices;
    for (size_t i = 0; i < n; i++) { rotIndices.push_back(i); rotIndices.push_back(-i); }
    cryptoContext->EvalRotateKeyGen(keyPair.secretKey, rotIndices);
    // Populate column containers with encryptions of isolated matrix elements.
    std::vector<std::vector<Ciphertext<DCRTPoly>>> enc_col_container; 
    for (size_t row=0 ; row < n ; ++row){ 
        for (size_t elem=0 ; elem < n ; ++elem){ 
            // Isolate row element and shift element to corresponding position in column.
            auto masked_enc_row = cryptoContext->EvalMult(enc_rows[row], encMasks[elem]); // Masked enc(row).
            cryptoContext->ModReduceInPlace(masked_enc_row);
            auto enc_elem = cryptoContext->EvalRotate(masked_enc_row, elem - row);
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
    std::vector<Ciphertext<DCRTPoly>> enc_cols; 
    for (size_t col=0 ; col < n ; ++col){ 
        enc_cols.push_back(cryptoContext->EvalAddMany(enc_col_container[col]));
    }   
    return enc_cols;
}

std::vector<Ciphertext<DCRTPoly>> evalMatrixSqMul(std::vector<Ciphertext<DCRTPoly>> &enc_rows, 
                                                int exponent,
                                                CryptoContext<DCRTPoly> &cryptoContext,
                                                InitMatrixExp initMatrixExp) {
    // TODO: init class for key rotation etc.
    // TODO: generate encrypted columns.
    // TODO: Square and multiply.
}


class InitMatrixExp {
public:
    InitMatrixExp(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair, int slots) :
        slots(slots) 
    {
        // TODO: assert slots leq plaintext slots in cryptoContext.
        // Generate rotation keys for |slots| steps.
        std::vector<int32_t> rotIndices;
        for (size_t i = 0; i <= slots; i++) { rotIndices.push_back(-i); }
        cryptoContext->EvalRotateKeyGen(keyPair.secretKey, rotIndices);
        // Generate encrypted masks.
        for (size_t elem=0 ; elem < slots ; ++elem){ 
            // Isolate row element and shift element to corresponding position in column.
            std::vector<int64_t> mask(slots,0); mask[elem] = 1; // Element mask.
            encMasks_.push_back(cryptoContext->Encrypt(keyPair.publicKey,
                                                   cryptoContext->MakePackedPlaintext(mask)));
        }
    }

    std::vector<Ciphertext<DCRTPoly>> encMasks() { return encMasks_; }

    const int slots;

private:
    std::vector<Ciphertext<DCRTPoly>> encMasks_;
};


// Matrix exponentiation of n x n matrix, encrypted by rows. Multiplicative depth: log(exponent)
// (Key required to generate rotation keys).
std::vector<Ciphertext<DCRTPoly>> evalMatrixExp(std::vector<Ciphertext<DCRTPoly>> &enc_rows, 
                                            int exponent,
                                            CryptoContext<DCRTPoly> &cryptoContext,
                                            KeyPair<DCRTPoly> keyPair,
                                            InitMatrixExp initMatrixExp) {

    auto matrix_dim = enc_rows.size();
    assert(initMatrixExp.slots == matrix_dim);

    // Generate enc(cols) - rowToColEnc generates rotation keys for context.
    auto enc_cols = rowToColEnc(enc_rows, cryptoContext, keyPair, initMatrixExp.encMasks());

    // Compute all encrypted rows and columns shifted by 0,1,...,n indices.
    std::vector<std::vector<Ciphertext<DCRTPoly>>> enc_rows_shifted;
    std::vector<std::vector<Ciphertext<DCRTPoly>>> enc_cols_shifted;
    for (size_t idx=0 ; idx < matrix_dim ; ++idx){ 
        auto enc_row_copy = cryptoContext->EvalRotate(enc_rows[idx], -matrix_dim); // enc(0...0|row)  
        auto enc_col_copy = cryptoContext->EvalRotate(enc_cols[idx], -matrix_dim); // enc(0...0|col)
        enc_rows[idx] = cryptoContext->EvalAdd(enc_rows[idx], enc_row_copy); // enc(row|row)  
        enc_cols[idx] = cryptoContext->EvalAdd(enc_cols[idx], enc_col_copy); // enc(col|col)  
        std::vector<Ciphertext<DCRTPoly>> enc_row_shifted;
        std::vector<Ciphertext<DCRTPoly>> enc_col_shifted;
        for (size_t shift=0 ; shift < matrix_dim ; ++shift){    
            enc_row_shifted.push_back(cryptoContext->EvalRotate(enc_rows[idx], shift));
            enc_col_shifted.push_back(cryptoContext->EvalRotate(enc_cols[idx], shift));
        }
        enc_rows_shifted.push_back(enc_row_shifted);
        enc_cols_shifted.push_back(enc_col_shifted);
    }
    // n x n Matrix exponentiation. Compute result row by row.
    std::vector<Ciphertext<DCRTPoly>> enc_rows_res; 
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
        int indices_dim = exponent-1; vectorIter indices(matrix_dim,indices_dim); bool iterContinue = true;  
        while (iterContinue) { 
            // std::cout << indices.value() << std::endl;
            std::vector<Ciphertext<DCRTPoly>> enc_mult_container; // container for factors of multiplicative terms.
            
            for (int i = 1; i < (indices_dim-2); ++i) { 
                // Push (i, j)  , (i+1, j),   (i+2, j) to container.   
                enc_mult_container.push_back(enc_cols_shifted[indices.value()[i]][indices.value()[i-1]]); 
                // Push (j, k) ,  (j, k+1)  , (j, k+2) to container.   
                enc_mult_container.push_back(enc_rows_shifted[indices.value()[i]][indices.value()[i+1]]);
            }
            // Push (k, l)  , (k+1, l)  , (k+2, l) to container.
            enc_mult_container.push_back(enc_cols_shifted[indices.value().back()][indices.value()[indices_dim-2]]); 
            // Push (l, 0) ,  (l, 1)    , (l, 2) to container.
            enc_mult_container.push_back(enc_rows_shifted[indices.value().back()][0]);
            // Multiply elements in enc_mult_container, add ciphertext to the enc_add_container.
            auto enc_mult = cryptoContext->EvalMultMany(enc_mult_container);

           // Multiply enc_mult with (row, i), (row, i+1), (row, i+2) ... for all rows.
            for (size_t row=0 ; row < matrix_dim ; ++row){
                enc_add_per_row_container[row].push_back(cryptoContext->EvalMult(enc_mult,
                                                                                 enc_rows_shifted[row][indices.value()[0]]));
            }

            iterContinue = indices.iterate(); // Continue iterating through i,j,k,l ...
        }

        // Sum all additive terms in enc_add_per_row_container[row], add result to each row encryption.
        for (size_t row=0 ; row < matrix_dim ; ++row){  
            enc_rows_res.push_back(cryptoContext->EvalAddMany(enc_add_per_row_container[row])); 
        }
        
        std::cout << "completed." << std::endl;
    // }
    return enc_rows_res;
}


class InitMatrixVecMult {
public:
    InitMatrixVecMult(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair, int slots) :
        slots(slots) 
    {
        // TODO: assert slots leq plaintext slots in cryptoContext.
        // Generate rotation keys for |slots| steps.
        std::vector<int32_t> rotIndices;
        for (size_t i = 0; i < slots; i++) { rotIndices.push_back(-i); }
        cryptoContext->EvalRotateKeyGen(keyPair.secretKey, rotIndices);
        // Generate Eval Sum Key for EvalInnerProduct.
        cryptoContext->EvalSumKeyGen(keyPair.secretKey);
        // Generate encrypted mask with leading 1.
        std::vector<int64_t> mask(slots,0); mask[0] = 1;
        encMaskFirst_ = cryptoContext->Encrypt(keyPair.publicKey, cryptoContext->MakePackedPlaintext(mask));
    }

    Ciphertext<DCRTPoly> encMaskFirst() { return encMaskFirst_; }

    const int slots;

private:
    Ciphertext<DCRTPoly> encMaskFirst_;
};


Ciphertext<DCRTPoly> evalMatrixVecMult(std::vector<Ciphertext<DCRTPoly>> &enc_rows, 
                                       Ciphertext<DCRTPoly> &enc_vec,
                                       CryptoContext<DCRTPoly> &cryptoContext,            
                                       InitMatrixVecMult &initMatrixVecMult) {
    // TODO: assert slots in initMatrixVecMult consistent with inputs.
    std::vector<Ciphertext<DCRTPoly>> enc_elements;
    for (size_t row=0 ; row < enc_rows.size() ; ++row){ 
        auto enc_element = cryptoContext->EvalInnerProduct(enc_rows[row], enc_vec, 
                                                           enc_rows.size());
        auto enc_element_masked = cryptoContext->EvalMult(enc_element, initMatrixVecMult.encMaskFirst());         
        cryptoContext->ModReduceInPlace(enc_element_masked);
        enc_element = cryptoContext->EvalRotate(enc_element_masked,-row);
        enc_elements.push_back(enc_element);
    }
    return cryptoContext->EvalAddMany(enc_elements);
}

// TODO: Requires initMatrixVecMult, initMatrixExp. Move these to separate context.
Ciphertext<DCRTPoly> evalVecMatrixMult(Ciphertext<DCRTPoly> &enc_vec,
                                       std::vector<Ciphertext<DCRTPoly>> &enc_rows, 
                                       KeyPair<DCRTPoly> keyPair, // TODO: Only public key required for row to col encryption.
                                       CryptoContext<DCRTPoly> &cryptoContext,            
                                       InitMatrixVecMult &initMatrixVecMult,
                                       InitMatrixExp &initMatrixExp) {
    // TODO: assert slots in initMatrixVecMult consistent with inputs.
    std::vector<Ciphertext<DCRTPoly>> enc_elements;
    for (size_t row=0 ; row < enc_rows.size() ; ++row){ 
        // Require col encryptions of matrix.
        auto enc_cols = rowToColEnc(enc_rows,cryptoContext, keyPair, initMatrixExp.encMasks()); 
        auto enc_element = cryptoContext->EvalInnerProduct(enc_cols[row], enc_vec, 
                                                           enc_cols.size());
        auto enc_element_masked = cryptoContext->EvalMult(enc_element, initMatrixVecMult.encMaskFirst());         
        cryptoContext->ModReduceInPlace(enc_element_masked);
        enc_element = cryptoContext->EvalRotate(enc_element_masked,-row);
        enc_elements.push_back(enc_element);
    }
    return cryptoContext->EvalAddMany(enc_elements);
}

// Ciphertext exponentiation, via square and multiply. Multiplicative depth: log(exponent)
Ciphertext<DCRTPoly> evalExponentiate(Ciphertext<DCRTPoly> &ciphertext, int exponent, 
                                      CryptoContext<DCRTPoly> &cryptoContext) {
    // Get msb position of p.
    int numBits = sizeof(int) * 8;
    int msbPosition = -1;
    for (int i = numBits - 1; i >= 0; i--) {
        if ((exponent >> i & 1) && (msbPosition == -1)) {
            msbPosition = i + 1;
            break;
        }
    }
    // Compute squarings of p.
    std::vector<Ciphertext<DCRTPoly>> ciphertexts_squarings;
    ciphertexts_squarings.push_back(ciphertext);
    for (int i = 1; i < msbPosition; i++) {
        ciphertexts_squarings.push_back(cryptoContext->EvalMult(ciphertexts_squarings[i-1], 
                                                                ciphertexts_squarings[i-1]));
    }
    // Select required squarings.
    std::vector<Ciphertext<DCRTPoly>> ciphertexts_squarings_container;
    for (int i = 0; i < msbPosition; i++) {
        if (exponent >> i & 1) {
            ciphertexts_squarings_container.push_back(ciphertexts_squarings[i]) ;
        }
    }
    // Multiply selected squarings.
    return cryptoContext->EvalMultMany(ciphertexts_squarings_container);
}


class InitPrefixMult {
public:
    InitPrefixMult(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair, int slots) :
        slots(slots) 
    {
        // TODO: assert slots leq plaintext slots in cryptoContext
        // Generate rotation keys.
        std::vector<int32_t> rotIndices;
        for (size_t i = 0; i < slots; i++) { rotIndices.push_back(-i); }
        cryptoContext->EvalRotateKeyGen(keyPair.secretKey, rotIndices);
    }

    const int slots;
};

Ciphertext<DCRTPoly> evalPrefixMult(Ciphertext<DCRTPoly> &ciphertext,
                                   int slots, CryptoContext<DCRTPoly> &cryptoContext,
                                   InitPrefixMult &initPrefixMult) {
    // Assert slots leq rotation keys. 
    assert(initPrefixMult.slots >= slots);

    // Prefix computation: 
    int levels = std::ceil(std::log2(slots));
        
    // Precompute rotation steps and plaintext masks. 
    // TODO: Move to precomputation (low priority, no cryptographic ops).
    std::vector<int32_t> rotSteps; std::vector<Plaintext> leadingOnesPlaintxts;
    for (size_t i = 0; i < levels; i++) { 
        rotSteps.push_back(std::pow(2, i)); 
        std::vector<int64_t> prefixOnes(slots,0); 
        for (size_t j = 0; j < rotSteps.back(); j++) { prefixOnes[j] = 1; }
        leadingOnesPlaintxts.push_back(cryptoContext->MakePackedPlaintext(prefixOnes));
    }

    // Compute prefix multiplications.
    auto ciphertext1 = ciphertext;
    for (size_t i = 0; i < levels; i++) {
        auto ciphertext2 = cryptoContext->EvalRotate(ciphertext1, -rotSteps[i]);
        // Pad ciphertext2 with leading 1's.        
        ciphertext2 = cryptoContext->EvalAdd(ciphertext2, leadingOnesPlaintxts[i]);
        // multiply ciphertext1 and ciphertext 2
        ciphertext1 = cryptoContext->EvalMult(ciphertext1, ciphertext2);
        cryptoContext->ModReduceInPlace(ciphertext1);
    }
    return ciphertext1;
} 

class InitPreserveLeadOne {
public:
    InitPreserveLeadOne(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair, int slots) :
        initPrefixMult(cryptoContext, keyPair, slots), slots(slots) 
    {
        // TODO: assert slots leq plaintext slots in cryptoContext
        // Generate rotation keys.
        cryptoContext->EvalRotateKeyGen(keyPair.secretKey, {-1});
        // Initialize prefix multiplication.
        // Generate encryption of masks. 
        std::vector<int64_t> ones(slots,1);
        std::vector<int64_t> leadingOne(slots,0); leadingOne[0] = 1;
        std::vector<int64_t> negOnes(slots,cryptoContext->GetCryptoParameters()->GetPlaintextModulus()-1);
        encOnes_ = cryptoContext->Encrypt(keyPair.publicKey,
                                          cryptoContext->MakePackedPlaintext(ones));
        encNegOnes_ = cryptoContext->Encrypt(keyPair.publicKey,
                                             cryptoContext->MakePackedPlaintext(negOnes));
        encLeadingOne_ = cryptoContext->Encrypt(keyPair.publicKey,
                                                cryptoContext->MakePackedPlaintext(leadingOne));
    }

    Ciphertext<DCRTPoly> encOnes() { return encOnes_; }
    Ciphertext<DCRTPoly> encNegOnes() { return encNegOnes_; }
    Ciphertext<DCRTPoly> encLeadingOne() { return encLeadingOne_; }

    InitPrefixMult initPrefixMult;
    const int slots;

private:
    Ciphertext<DCRTPoly> encOnes_;
    Ciphertext<DCRTPoly> encNegOnes_;
    Ciphertext<DCRTPoly> encLeadingOne_;
};

Ciphertext<DCRTPoly> evalPreserveLeadOne(Ciphertext<DCRTPoly> &ciphertext,
                                         CryptoContext<DCRTPoly> &cryptoContext,
                                         InitPreserveLeadOne &initPreserveLeadOne) {
    // (1-x0),(1-x1),...,(1-xn).
    auto encDiffs = cryptoContext->EvalAdd(cryptoContext->EvalMult(ciphertext,initPreserveLeadOne.encNegOnes()),
                                           initPreserveLeadOne.encOnes());
    // y0, y1,..., yn: yi = ith multiplicative prefix.
    auto encPrefix = evalPrefixMult(encDiffs, initPreserveLeadOne.slots, 
                                    cryptoContext, initPreserveLeadOne.initPrefixMult);
    // x0, x1*y0 ,...,   xn*yn-1
    auto result = cryptoContext->EvalMult(ciphertext,
                                   cryptoContext->EvalAdd(initPreserveLeadOne.encLeadingOne(), 
                                   cryptoContext->EvalRotate(encPrefix,-1)));
    cryptoContext->ModReduceInPlace(result);
    return result;
}

class InitNotEqualZero {
public:
    InitNotEqualZero(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair, int slots, int range ) :
        slots(slots), range(range)
    {
        // TODO: assert range leq plaintext slots in cryptoContext
        // Precompute constants.
        auto plaintxtModulus = cryptoContext->GetCryptoParameters()->GetPlaintextModulus();
        int factorialRange = modFactorial(range, plaintxtModulus); 
        int invFactorialRange = modInverse(factorialRange, plaintxtModulus); 

        // Packed encryption of constants.
        std::vector<int64_t> invfPacked(slots,invFactorialRange);
        encInvFactorial_= cryptoContext->Encrypt(keyPair.publicKey,
                                                 cryptoContext->MakePackedPlaintext(invfPacked));
        std::vector<int64_t> onePacked(slots,1); 
        encOne_ = cryptoContext->Encrypt(keyPair.publicKey,
                                         cryptoContext->MakePackedPlaintext(onePacked));    
        for (size_t i=1 ; i <= range ; ++i){ 
        std::vector<int64_t> negIntPacked(slots,plaintxtModulus-i); 
        encNegRange_.push_back(cryptoContext->Encrypt(keyPair.publicKey,
                              cryptoContext->MakePackedPlaintext(negIntPacked)));
        };
    }

    Ciphertext<DCRTPoly> encOne() { return encOne_; }
    Ciphertext<DCRTPoly> encInvFactorial() { return encInvFactorial_; }
    std::vector<Ciphertext<DCRTPoly>> encNegRange() { return encNegRange_; }

    const int slots;
    const int range;

private:
    Ciphertext<DCRTPoly> encOne_;
    Ciphertext<DCRTPoly> encInvFactorial_;
    std::vector<Ciphertext<DCRTPoly>> encNegRange_;
};

Ciphertext<DCRTPoly> evalNotEqualZero(Ciphertext<DCRTPoly> &ciphertext,
                                  CryptoContext<DCRTPoly> &cryptoContext,
                                  InitNotEqualZero &initNotEqualZero) {
    // Compute binary map for range r: 1-(x-1)(x-2)...(x-r)/r!
    std::vector<Ciphertext<DCRTPoly>> encDiffs;
    for (size_t i=0 ; i < initNotEqualZero.range ; ++i){ 
        encDiffs.push_back(cryptoContext->EvalAdd(ciphertext, initNotEqualZero.encNegRange()[i]));
    }
    encDiffs.push_back(initNotEqualZero.encInvFactorial());
    auto encMult = cryptoContext->EvalMultMany(encDiffs);
    if (initNotEqualZero.range % 2) { return cryptoContext->EvalAdd(initNotEqualZero.encOne(),encMult); } 
    else { return cryptoContext->EvalAdd(initNotEqualZero.encNegRange()[0],encMult); }
}

void refreshInPlace(Ciphertext<DCRTPoly> &ciphertext, int slots, 
                    KeyPair<DCRTPoly> keyPair, CryptoContext<DCRTPoly> &cryptoContext){
    Plaintext plaintextExpRes;
    cryptoContext->Decrypt(keyPair.secretKey, ciphertext, &plaintextExpRes); 
    plaintextExpRes->SetLength(slots); auto payload = plaintextExpRes->GetPackedValue();
    ciphertext = cryptoContext->Encrypt(keyPair.publicKey, cryptoContext->MakePackedPlaintext(payload));
}



int main(int argc, char* argv[]) {

    ////////////////////////////////////////////////////////////
    // Set-up of parameters
    ////////////////////////////////////////////////////////////

    // benchmarking variables
    TimeVar t;
    double processingTime(0.0);

    // Crypto Parameters
    // # of evalMults = 3 (first 3) is used to support the multiplication of 7
    // ciphertexts, i.e., ceiling{log2{7}} Max depth is set to 3 (second 3) to
    // generate homomorphic evaluation multiplication keys for s^2 and s^3
    CCParams<CryptoContextBGVRNS> parameters;

    int chosen_ptxtmodulus = 65537;
    parameters.SetPlaintextModulus(chosen_ptxtmodulus);
    // p = 65537, depth = 13 -> "Please provide a q and a m satisfying: (q-1)/m is an integer. The values of primeModulus = 65537 and m = 131072 do not."
    // Fermats thm works for p = 786433, dep = 20.
    parameters.SetMultiplicativeDepth(12);
    parameters.SetMaxRelinSkDeg(3);


    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    std::cout << "Ring dimension N: " << cryptoContext->GetRingDimension() << std::endl;

    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);

    std::cout << "Plaintext modulus p = " << cryptoContext->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
    std::cout << "Cyclotomic order n = " << cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2
              << std::endl;
    // std::cout << "log2 q = "
    //           << log2(cryptoContext->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
    //           << std::endl;

    // Initialize Public Key Containers
    KeyPair<DCRTPoly> keyPair;

    // Perform Key Generation Operation
    TIC(t);

    keyPair = cryptoContext->KeyGen();

    processingTime = TOC(t);
    std::cout << "Key generation time: " << processingTime << "ms" << std::endl;

    if (!keyPair.good()) {
        std::cout << "Key generation failed!" << std::endl;
        exit(1);
    }

    std::cout << "Running key generation for homomorphic multiplication "
                 "evaluation keys..."
              << std::endl;

    TIC(t);

    cryptoContext->EvalMultKeysGen(keyPair.secretKey);

    processingTime = TOC(t);
    std::cout << "Key generation time for homomorphic multiplication evaluation keys: " << processingTime << "ms"
              << std::endl;


    ////////////////////////////////////////////////////////////
    // Top Trading Cycle Algorithm.
    ////////////////////////////////////////////////////////////

    //==========================================================
    // Offline phase.
    //==========================================================

    // Offline: User preferences.
    std::vector<std::vector<int64_t>> userInputs;
    userInputs.push_back({4, 1, 2, 3, 0});
    userInputs.push_back({4, 3, 2, 1, 0});
    userInputs.push_back({4, 1, 0, 2, 3});
    userInputs.push_back({1, 3, 4, 0, 2});
    userInputs.push_back({3, 1, 2, 0, 4});
    auto n = userInputs.size();

    // Offline: User preference as a permutation matrix.
    std::vector<std::vector<Ciphertext<DCRTPoly>>> encUserPrefList;
    std::vector<std::vector<Ciphertext<DCRTPoly>>> encUserPrefTransposedList;
    for (int user=0; user < n; ++user) { 
        std::vector<std::vector<int64_t>> userPrefMatrix;
        for (int j=0; j < n; ++j) {
            std::vector<int64_t> row(n,0); row[userInputs[user][j]] = 1;
            userPrefMatrix.push_back(row);
        }   
        // Transpose user preference-permutation matrix.
        std::vector<std::vector<int64_t>> userPrefMatrixTransposed;
        std::vector<int64_t> zeroRow(n,0);
        for (int j=0; j < n; ++j) { userPrefMatrixTransposed.push_back(zeroRow); }
        for (int j=0; j < n; ++j) { for (int k=0; k < n; ++k) {
            if (userPrefMatrix[j][k] == 1) { userPrefMatrixTransposed[k][j] = 1; break; }
        }}
        // Encrypt user preference and transposed preference permutation matrix.
        std::vector<Ciphertext<DCRTPoly>> encUserPref;
        std::vector<Ciphertext<DCRTPoly>> encUserPrefTransposed;
        for (int j=0; j<n ; ++j){ 
            encUserPref.push_back(cryptoContext->Encrypt(keyPair.publicKey,
                                  cryptoContext->MakePackedPlaintext(userPrefMatrix[j])));
            encUserPrefTransposed.push_back(cryptoContext->Encrypt(keyPair.publicKey,
                                            cryptoContext->MakePackedPlaintext(userPrefMatrixTransposed[j])));
        }
        encUserPrefList.push_back(encUserPref);
        encUserPrefTransposedList.push_back(encUserPrefTransposed);
    }

    // Server Offline: Precompute encryptions of constants, initialize rotation keys.
    InitMatrixExp initMatrixExp(cryptoContext, keyPair, userInputs.size());
    InitMatrixVecMult initMatrixVecMult(cryptoContext, keyPair, userInputs.size());
    InitPrefixMult initPrefixMult(cryptoContext, keyPair, userInputs.size());
    InitNotEqualZero initNotEqualZero(cryptoContext, keyPair, userInputs.size(), userInputs.size()); // (..., slots, range)
    InitPreserveLeadOne initPreserveLeadOne(cryptoContext,keyPair,n);
    
    // Server Offline: Initialize enc(user-availability), enc(ones), enc([0:n])
    // TODO: move to initialization class.
    std::vector<int64_t> ones(n,1);
    std::vector<int64_t> negOnes(n,cryptoContext->GetCryptoParameters()->GetPlaintextModulus()-1);
    std::vector<int64_t> range; for (int i=0; i<n; ++i) { range.push_back(i); }
    auto encUserAvailability = cryptoContext->Encrypt(keyPair.publicKey,
                                                      cryptoContext->MakePackedPlaintext(ones));
    auto encOnes = encUserAvailability; 
    auto encNegOnes = cryptoContext->Encrypt(keyPair.publicKey,
                                             cryptoContext->MakePackedPlaintext(negOnes));
    auto encRange = cryptoContext->Encrypt(keyPair.publicKey,
                                           cryptoContext->MakePackedPlaintext(range));
                                                  
    //==========================================================
    // Online phase.
    //==========================================================

    // Init output ciphertext.
    auto enc_output = encNegOnes; // -1 output => not on cycle.

    // TODO: loop over n rounds.

    //----------------------------------------------------------
    // (1) Update adjacency matix.
    //----------------------------------------------------------
    
    TIC(t);
    auto res = cryptoContext->EvalMult(encOnes,encOnes);
    processingTime = TOC(t);
    std::cout << "Multiplication: " << processingTime << "ms" << std::endl;

    TIC(t);

    std::vector<Ciphertext<DCRTPoly>> encRowsAdjMatrix;

    // Generate adjacency matrix for all users.
    for (int user = 0; user < n; ++user){
        // Sort availability according to user preference.
        auto encUserAvailablePref = evalMatrixVecMult(encUserPrefList[user], encUserAvailability,
                                                      cryptoContext, initMatrixVecMult);
        // Preserve highest, available preference.
        auto encUserFirstAvailablePref = evalPreserveLeadOne(encUserAvailablePref, cryptoContext, initPreserveLeadOne);

        // Transpose back to obtain adjacency matrix row.
        encRowsAdjMatrix.push_back(evalMatrixVecMult(encUserPrefTransposedList[user], encUserFirstAvailablePref,
                                                       cryptoContext, initMatrixVecMult));
    }

    std::cout << "Adjacency Matrix: " << std::endl;
    for (int row=0; row < n; ++row){
        Plaintext plaintext;
        cryptoContext->Decrypt(keyPair.secretKey, encRowsAdjMatrix[row], &plaintext); 
        plaintext->SetLength(n); auto payload = plaintext->GetPackedValue();
        std::cout << payload << std::endl;
    }

    processingTime = TOC(t);
    std::cout << "Online part 1 - Adjacency matrix update time: " << processingTime << "ms" << std::endl;

    // Refresh ciphertexts.
    for (int row=0; row < n; ++row){ refreshInPlace(encRowsAdjMatrix[row],n, keyPair, cryptoContext); } 

    //----------------------------------------------------------
    // (2) Matrix exponentiation for cycle finding.
    //----------------------------------------------------------

    TIC(t);

    // TODO: get compute encryption of columns.
    // TODO: compute exponentiation (enc(rows),enc(cols))


    auto encMatrixExp = evalMatrixExp(encRowsAdjMatrix,n,cryptoContext,keyPair,initMatrixExp); // TODO: access keypair from InitMatrixExp class.
    
    // Refresh ciphertexts.
    for (int row=0; row < n; ++row){ refreshInPlace(encMatrixExp[row],n,keyPair, cryptoContext); } 

    // u: Derive computed cycle.
    auto enc_u = evalVecMatrixMult(encOnes,encMatrixExp,keyPair,cryptoContext,initMatrixVecMult,initMatrixExp); // TODO: only publickey required, create dedicated init class.
    enc_u = evalNotEqualZero(enc_u,cryptoContext,initNotEqualZero); 
    
    processingTime = TOC(t);
    std::cout << "Online part 2 - Matrix exponentiation time: " << processingTime << "ms" << std::endl;

    //----------------------------------------------------------
    // (3) Update user availability and outputs.
    //----------------------------------------------------------

    // Refresh ciphertexts.
    refreshInPlace(enc_u,n,keyPair, cryptoContext);

    TIC(t);

    // t: Compute current preference index for all users in packed ciphertext.
    std::vector<Ciphertext<DCRTPoly>> enc_elements;
    for (int user=0; user < n; ++user){ 
        auto enc_t_user = cryptoContext->EvalInnerProduct(encRowsAdjMatrix[user], encRange, 
                                                          encRowsAdjMatrix.size());
        enc_t_user = cryptoContext->EvalMult(enc_t_user, initMatrixVecMult.encMaskFirst()); // TODO: 
        cryptoContext->ModReduceInPlace(enc_t_user);
        enc_elements.push_back(cryptoContext->EvalRotate(enc_t_user,-user));
    }
    auto enc_t = cryptoContext->EvalAddMany(enc_elements);
    // o: Update output for all users in packed ciphertext: o <- t x u + o x (1-u)
    auto enc_t_mult_u = cryptoContext->EvalMult(enc_t, enc_u); cryptoContext->ModReduceInPlace(enc_t);
    auto enc_one_min_u = cryptoContext->EvalAdd(encOnes,
                                                cryptoContext->EvalMult(enc_u, encNegOnes));
    // output <- t x u + o x (1-u)
    enc_output = cryptoContext->EvalAdd(enc_t_mult_u, 
                                        cryptoContext->EvalMult(enc_output,enc_one_min_u));
    // Update availability: 1-NotEqualZero(output)
    auto enc_output_reduced = evalNotEqualZero(enc_output,cryptoContext,initNotEqualZero); 
    encUserAvailability = cryptoContext->EvalAdd(encOnes,
                                                 cryptoContext->EvalMult(enc_output_reduced, encNegOnes));

    processingTime = TOC(t);
    std::cout << "Online part 3 - User availability & output update time: " << processingTime << "ms" << std::endl;
    
    // Print output vector (-1 means not on cycle).
    Plaintext plaintext;
    cryptoContext->Decrypt(keyPair.secretKey, enc_output, &plaintext); 
    plaintext->SetLength(n); auto payload = plaintext->GetPackedValue();
    std::cout << "Cycle finding result: " << payload << std::endl;

    // TODO: end loop.
    
    return 0;
}
