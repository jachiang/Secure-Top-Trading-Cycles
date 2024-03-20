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

#include "utilities.h"
#include "crypto_utilities.h"
#include "crypto_enc_transform.h"
#include "crypto_matrix_operations.h"
#include "crypto_prefix_mult.h"
#include "crypto_noteqzero.h"

#include <cassert>
#include <iostream>
#include <iterator>
#include <cmath>
#include <vector>

#include "openfhe.h"

using namespace lbcrypto;


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
    InitRotsMasks initRotsMasks(cryptoContext,keyPair,n);
    // InitMatrixExp initMatrixExp(cryptoContext, keyPair, userInputs.size());
    // InitMatrixVecMult initMatrixVecMult(cryptoContext, keyPair, userInputs.size());
    // InitPrefixMult initPrefixMult(cryptoContext, keyPair, userInputs.size());
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
    // Tests for matrix exponentiation.
    //==========================================================
               
    // std::vector<std::vector<int64_t>> testMatRows;
    // testMatRows.push_back({1, 0, 1, 0, 0});
    // testMatRows.push_back({0, 1, 0, 0, 1});
    // testMatRows.push_back({1, 0, 1, 0, 0});
    // testMatRows.push_back({0, 1, 0, 0, 1});
    // testMatRows.push_back({1, 0, 1, 0, 0});

    // std::vector<Ciphertext<DCRTPoly>> testEncMatRows;
    // for (int i = 0; i < testMatRows.size(); i++){
    // testEncMatRows.push_back(cryptoContext->Encrypt(keyPair.publicKey,
    //                          cryptoContext->MakePackedPlaintext(testMatRows[i])));
    // }
    // Test 1:
    // auto encMatsResElems = evalMatSquarings(testEncMatRows,3,cryptoContext,initRotsMasks, keyPair);
    // auto encMatResElems = evalMatrixMul2Pow(encMatsResElems,cryptoContext,initRotsMasks);
    // auto encMatResRows = encElem2Rows(encMatResElems,cryptoContext,initRotsMasks);

    // Test 2:
    // auto encMatResElems = evalMatSqMul(testEncMatRows, 15, cryptoContext,initRotsMasks, keyPair);
    // auto encMatResRows = encElem2Rows(encMatResElems,cryptoContext,initRotsMasks);

    // std::cout << "Test Matrix: " << std::endl;
    // printEncMatRows(encMatResRows,cryptoContext,keyPair);

    // return 0;

    //==========================================================
    // Test for single multiplication.
    //==========================================================

    // TIC(t);
    // auto res = cryptoContext->EvalMult(encOnes,encOnes);
    // cryptoContext->ModReduceInPlace(res);
    // processingTime = TOC(t);
    // std::cout << "Multiplication & mod reduction: " << processingTime << "ms" << std::endl;
    // return 0;

    //==========================================================
    // Online phase.
    //==========================================================

    // Init output ciphertext.
    auto enc_output = encNegOnes; // -1 output => not on cycle.

    // TODO: loop over n rounds.

    //----------------------------------------------------------
    // (1) Update adjacency matix.
    //----------------------------------------------------------
    
    CryptoOpsLogger cryptoOpsLogger;


    TIC(t);

    std::vector<Ciphertext<DCRTPoly>> encRowsAdjMatrix;

    // Generate adjacency matrix for all users.
    for (int user = 0; user < n; ++user){
        // Sort availability according to user preference.
        auto encUserAvailablePref = evalMatrixVecMult(encUserPrefList[user], encUserAvailability,
                                                      cryptoContext, initRotsMasks);
        // Preserve highest, available preference.
        auto encUserFirstAvailablePref = evalPreserveLeadOne(encUserAvailablePref, cryptoContext, initPreserveLeadOne);

        // Transpose back to obtain adjacency matrix row.
        encRowsAdjMatrix.push_back(evalMatrixVecMult(encUserPrefTransposedList[user], encUserFirstAvailablePref,
                                                       cryptoContext, initRotsMasks));
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
    // auto encMatrixExp = evalMatrixExp(encRowsAdjMatrix,n,cryptoContext,initRotsMasks,cryptoOpsLogger); 
    auto encMatrixExpElems = evalMatSqMul(encRowsAdjMatrix,n,cryptoContext,initRotsMasks,cryptoOpsLogger,keyPair);
    auto encMatrixExp = encElem2Rows(encMatrixExpElems,cryptoContext,initRotsMasks,cryptoOpsLogger);
    std::cout << "Total matrix exponentiation time: " << TOC(t) << " ms" << std::endl;
    // std::cout << "Total homomorphic computation time: " << cryptoOpsLogger.totalTime() << " ms" << std::endl;
    // std::cout << "Total homomorphic multiplications: " << cryptoOpsLogger.multOps()+cryptoOpsLogger.innerProdOps() << std::endl;
    // std::cout << "Total homomorphic multiplication time: " << cryptoOpsLogger.multTime()+cryptoOpsLogger.innerProdTime() << " ms" << std::endl;

    // Refresh ciphertexts.
    for (int row=0; row < n; ++row){ refreshInPlace(encMatrixExp[row],n,keyPair, cryptoContext); } 

    // u: Derive computed cycle.
    auto enc_u = evalVecMatrixMult(encOnes,encMatrixExp,cryptoContext,initRotsMasks,cryptoOpsLogger); // TODO: only publickey required, create dedicated init class.
    enc_u = evalNotEqualZero(enc_u,cryptoContext,initNotEqualZero); 
    
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
        enc_t_user = cryptoContext->EvalMult(enc_t_user, initRotsMasks.encMasks()[0]); // TODO: 
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
