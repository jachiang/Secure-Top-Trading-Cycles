from functools import reduce

# Define matrix
m = [[1,2,3],
     [4,5,6],
     [7,8,9]]

# Memoization table
# M_{r,i} x M_{i,j} x M_{j,c}

# Determine sizes of elements
# 2^2 + 2^1 + 2^0

# Key for tables: 
# List of indices (i,j),(j,k),(k,l),(l,m) -> M_{i,j} M_{j,k} M_{k,l} M_{l,m} 
# To produce key, must sort list of indices (order doesn't matter).
# key <-- hash(tuple(sorted(indexList)))

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
            print("hit")
            return value 
        else:
            elemContainer = []
            # If number of elements is a power of two.
            if sum(binList) == 1:
                # If number of elements is 2.
                if len(indexPairList) == 2:
                    # Retrieve elements in matrix
                    leftElement = self.matrix[indexPairList[0][0]][indexPairList[0][1]]
                    rightElement = self.matrix[indexPairList[1][0]][indexPairList[1][1]]
                    res = leftElement * rightElement
                    newKey = hash(tuple(sorted(indexPairList)))
                    self.table[newKey] = res
                    return res
                # If number of elements is power of 2.
                else:
                    # recursive calls.
            # If number of elements is not a power of two.
            else: 
                # Bitscan over length of indexPairList
                indexPairListCopy = indexPairList
                for i in range(len(binList)):
                    if int(binStr[binListLen-1-i]):
                        # For active bit a non-zero index: 
                        if i == 0:
                            # Append single element to container.
                            elementIdxPair = indexPairListCopy[-1]
                            indexPairListCopy = indexPairListCopy[:-1]
                            elemContainer.append(self.matrix[elementIdxPair[0]][elementIdxPair[1]])
                        else: 
                            # Should break down 
                            elementIdxPairs = indexPairListCopy[-2**i:]
                            print(i, elementIdxPairs)
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

print(Test.retrieve([(1,1),(2,2),(1,2),(2,1)]))




# Input:
# ...