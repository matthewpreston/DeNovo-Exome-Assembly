def MergeLocations(locations):
    """Merges contiguous locations"""
    return MergeSortLocs(locations)
        
def MergeSortLocs(loc):
    """
    Uses merge sort algorithm to compress locations FULLY (I'm tired of using a
    linear sliding window algorithm only to have it fail repeatedly)
    """
    if len(loc) <= 1:
        return loc
    mid = len(loc) // 2
    left = loc[:mid]
    right = loc[mid:]
    
    left = MergeSortLocs(left)
    right = MergeSortLocs(right)
    
    result = []
    i = j = 0
    while i < len(left) and j < len(right):
        if left[i][1] <= right[j][1]:
            result.append(left[i])
            i += 1
        else:
            result.append(right[j])
            j += 1
 
    if left:
        result.extend(left[i:])
    if right:
        result.extend(right[j:])
    #print("Before compress:", result)
    result = CompressLocsFwd(result)
    #print("After compress:", result)
    return result

def CompressLocsFwd(loc, start=0):
    i = start
    while i+1 < len(loc):#Length may decrease
        #Test if endpoint of 1st tuple is within the range of the second
        #I.E. (1,4)+(5,6)->(1,6) since 5-1 <= 4 <= 6
        #print("Testing fwd", loc[i], loc[i+1])
        #print(loc[i+1][0]-1 <= loc[i][1])
        #print(loc[i][1] <= loc[i+1][1])
        if loc[i+1][0]-1 <= loc[i][1] and loc[i][1] <= loc[i+1][1]:
            temp = loc[i]
            loc[i] = (min(loc[i][0], loc[i+1][0]), loc[i+1][1])
            del loc[i+1]
            #print("During compress:", loc)
            #print("Testing temp", temp, loc[i])
            if temp[0] != loc[i][0]: #Changed, put it into reverse gear
                #print("Going backwards; i:", i)
                loc, i = CompressLocsBwd(loc, i)
                #print("Done backwards; i:", i)
        else:
            i += 1
    return loc

def CompressLocsBwd(loc, start):
    i = start
    while 0 <= i-1:
        #Test if endpoint of 1st tuple is within the range of the second
        #I.E. (1,4)+(5,6)->(1,6) since 5-1 <= 4 <= 6
        #print("Testing bwd", loc[i-1], loc[i])
        #print(loc[i][0]-1 <= loc[i-1][1])
        #print(loc[i-1][1] <= loc[i][1])
        if loc[i][0]-1 <= loc[i-1][1] and loc[i-1][1] <= loc[i][1]:
            temp = loc[i]
            loc[i] = (min(loc[i-1][0], loc[i][0]), loc[i][1])
            del loc[i-1]
            #print("During compress:", loc)
        else:
            break
        i -= 1
    return (loc, i)