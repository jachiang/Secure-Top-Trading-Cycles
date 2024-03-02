from functools import reduce

class MultLookup:
    def __init__(self, matrix):
        self.table = {}
        self.matrix = matrix
        self.computedMults = 0

    def resetCtr(self):
        self.computedMults = 0

    # Note: indexPairList is a list of int pairs
    # which reference position of a matrix element.
    def retrieveOrCompute(self, indexPairList):
        # Bitdecompose the number of pairs in indexPairList.
        binStr = bin(len(indexPairList))[2:] 
        binList = list(map(int, binStr))
        binListLen = len(binList)
        # Check if in memoization table.
        key = hash(tuple(sorted(indexPairList)))
        value = self.table.get(key)
        if value:
            self.tableHit = True
            return value 
        else:
            # If there is no hit, populate container with elements for multiplication:
            # Container elements are multiplicative terms of "powers-of-two" elements 
            # indicated in indexPairList. The multiplication of these terms in the container
            # is performed at the end of this branch.
            elemContainer = [] 
            # Case 0: If number of matrix elements referenced in indexPairList is a power of two.
            # e.g M_i,j M_j,k M_k,l M_l,m (4 matrix elements = 2^2)
            if sum(binList) == 1:
                # Case 0a: If number of matrix elements is 2.
                if len(indexPairList) == 2:
                    # Retrieve elements in matrix
                    leftElement = self.matrix[indexPairList[0][0]][indexPairList[0][1]]
                    rightElement = self.matrix[indexPairList[1][0]][indexPairList[1][1]]
                # Case 0b: If number of matrix elements is (greater than 1) power of 2.
                else:
                    indexPairListLeft = indexPairList[:int(len(indexPairList)/2)]
                    indexPairListRight = indexPairList[int(len(indexPairList)/2):]
                    leftElement = self.retrieveOrCompute(indexPairListLeft)
                    rightElement =  self.retrieveOrCompute(indexPairListRight)
                elemContainer.append(leftElement)
                elemContainer.append(rightElement)
            # Case 1: If number of matrix elements referenced in indexPairList is *not* power of two.
            # e.g M_i,j M_j,k M_k,l M_l,m M_m,n (5 matrix elements = 2^2+2^0)
            else: 
                # Split referenced matrix elements into batches of 2-power sizes.
                indexPairListCopy = indexPairList
                for i in range(len(binList)):
                    if int(binStr[binListLen-1-i]):
                        if i == 0:
                            # Append single element to container.
                            elementIdxPair = indexPairListCopy[-1]
                            indexPairListCopy = indexPairListCopy[:-1]
                            elemContainer.append(self.matrix[elementIdxPair[0]][elementIdxPair[1]])
                        else: 
                            # Append batch of elements (of the 2-power size) to container.
                            elementIdxPairs = indexPairListCopy[-2**i:]
                            indexPairListCopy = indexPairListCopy[:-2**i]
                            elemContainer.append(self.retrieveOrCompute(elementIdxPairs)) 
            # Multiply elements in elemContainer
            self.computedMults += len(elemContainer)-1
            res = reduce(lambda x, y: x * y, elemContainer)
            # Add to hash table.
            newKey = hash(tuple(sorted(indexPairList)))
            self.table[newKey] = res
            return res

class indicesIter:
    def __init__(self, modulus, slots):
        self.mod = modulus        
        self.indices = []
        for i in range(slots):
            self.indices.append(0)
        self.computedMults = 0

    def step(self):
        return self.step_(0)
    
    def step_(self, slot):
        # Case 0: Highest idx cannot be incremented further.
        if slot == len(self.indices)-1 and self.indices[slot] == self.mod-1:
            return False
        # Case 1: Increment & carry. 
        elif self.indices[slot] == self.mod -1:
            self.indices[slot] = 0
            if self.step_(slot+1):
                return True
            else:
                return False
        # Case 2: Increment.
        else:
            self.indices[slot] =  self.indices[slot]+1
            return True



def main():

    # Set matrix to be exponentiated.
    matrix = [[1,2,3],
              [4,5,6],
              [7,8,9]]
    exp = 5

    # Multiplications computed without lookup.
    multsNoLookup = 0

    # Initialize resulting matrix.
    res = []
    for i in range(len(matrix)):
        col = []
        for i in range(len(matrix)):
            col.append(0)
        res.append(col)

    # Initialize lookup table for multiplications
    multLookup = MultLookup(matrix)

    # Compute matrix exponentiation element-wise.
    matrixIter = indicesIter(len(matrix),2)
    contMatrixIter = True
    while(contMatrixIter):
        (r,c) = matrixIter.indices[0], matrixIter.indices[1]
        # Fill container with mult. terms for summation over
        # all i,j,k,l in sum  M_ri M_ij M_jk M_kl M_lc 
        sumContainer = []
        sumIter = indicesIter(len(matrix),exp-1)
        contSummation = True
        while(contSummation):
            # Compute multiplication M_ri M_ij M_jk M_kl M_lc 
            # Construct index pairs (r,i), (i,j), ...
            indexPairList = []
            indexPairList.append((r,sumIter.indices[0]))
            for slot in range(len(sumIter.indices)-1):
                idxPair = (sumIter.indices[slot],sumIter.indices[slot+1])
                indexPairList.append(idxPair)
            indexPairList.append((sumIter.indices[-1],c))
            # Add to sumContainer
            sumContainer.append(multLookup.retrieveOrCompute(indexPairList))
            # Capture mults required
            multsNoLookup += exp-1
            contSummation = sumIter.step()
        # Compute final summation for matrix term M_rc
        res[r][c] = reduce(lambda x, y: x + y, sumContainer)
        contMatrixIter = matrixIter.step()

    print(res)
    print("Mults with lookup:", multLookup.computedMults)
    print("Mults without lookup:", multsNoLookup)


if __name__ == "__main__":
    main()