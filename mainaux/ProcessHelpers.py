import operator
import numpy as np

def roll_dice(num, rate):
    # if rate is > 1
    # switch to Poisson
    if rate > 1:
        temp_num = np.random.poisson(num * rate)
        rtn_amt = min(temp_num, num)
    else:
        # rate is <= 1
        temp_rand = np.random.rand(num) # creates an array of rand elem of size 1 x num
        rtn_amt = (temp_rand < rate).sum()
    return rtn_amt

# moves amt from input_dict[key1] to input_dict[key2]
# Assumes amt is a valid amount
def move_buckets(input_dict, key1, key2, amt):
    input_dict[key1] -= amt
    input_dict[key2] += amt

def transfer_buckets(input_dict, key1, key2, rate):
    amt_to_transfer = roll_dice(input_dict[key1], rate)
    move_buckets(input_dict, key1, key2, amt_to_transfer)
    return amt_to_transfer

ops =   {
            "<" : operator.lt,
            "<=": operator.le,
            "==": operator.eq,
            "!=": operator.ne,
            ">=": operator.ge,
            ">" : operator.gt
        }    