#Uses the current installed version of RAT to build the current data
#Cleaning bitmask.  Can also
import numpy as np
import rat

def get_dcwords():
    dcbits = rat.RAT.DU.DataCleaningBits()
    dcbits.BeginOfRun()
    #Hack to find all bits
    dcdict = {}
    bitindex = 0
    maxbit = 63
    for bit in range(maxbit):
        try:
            name = dcbits.GetBitName(bit)
            dcdict[bit] = name
        except:
            break
    return dcdict

def binary_bit_to_int(bitindex):
    return 2**bitindex

def int_to_binary_bit(integer):
    return np.log(integer)/np.log(2)

if __name__=='__main__':
    test_dict = get_dcwords()
    print(test_dict)
