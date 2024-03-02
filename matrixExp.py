from functools import reduce

class MultLookup:
    def __init__(self, matrix):
        self.table = {}
        self.matrix = matrix

    # Note: indexPairList is a list of int pairs.
    def retrieve(self, indexPairList):
        # Bitdecompose the number of pairs in indexPairList.
        binStr = bin(len(indexPairList))[2:] 
        binList = list(map(int, binStr))
        binListLen = len(binList)
        # Check if in memoization table.
        key = hash(tuple(sorted(indexPairList)))
        value = self.table.get(key)
        if value:
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
                    leftElement = self.retrieve(indexPairListLeft)
                    rightElement =  self.retrieve(indexPairListRight)
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
                            elemContainer.append(self.retrieve(elementIdxPairs)) 
            # Multiply elements in elemContainer
            res = reduce(lambda x, y: x * y, elemContainer)
            # Add to hash table.
            newKey = hash(tuple(sorted(indexPairList)))
            self.table[newKey] = res
            return res
            
Test = MultLookup([[1,2,3],
                   [4,5,6],
                   [7,8,9]])

# 5, 9, 6, 8
print(Test.retrieve([(1,1),(2,2),(1,2),(2,1),(0,0)]))




# Input:
# ...