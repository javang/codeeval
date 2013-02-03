#!/usr/bin/python

import sys
import os
import logging
log = logging.getLogger("assembler")
from math import log10

ones = {  1: "One",
          2: "Two",
          3: "Three",
          4: "Four",
          5: "Five",
          6: "Six",
          7: "Seven",
          8: "Eight",
          9: "Nine",
         }
teens  = { 10: "Ten",
           11: "Eleven",
           12: "Twelve",
           13: "Thirteen",
           14: "Fourteen",
           15: "Fiveteen",
           16: "Sixteen",
           17: "Seventeen",
           18: "Eighteen",
           19: "Nineteen",
        }

tens = { 2: "Twenty",
           3: "Thirty",
           4: "Forty",
           5: "Fifty",
           6: "Sixty",
           7: "Seventy",
           8: "Eighty",
           9: "Ninety",
        }

powers = {0:"", 1: "Thousand", 2:"Million",3: "Billion"}

def number_to_text(number):
    if number == 0:
        return "ZeroDollars"
    words = []
    txt = str(number)
    i = 0
    while len(txt) > 0:
        n = int(txt[-3:])
        if n != 0:
            words.append( get3digits(n) + powers[i] )
        txt = txt[:-3]
        i+=1
    words.reverse()
    return "".join(words) + "Dollars"

def get1digit(val):
    if val > 0:
        return ones[val]
    return ""

def get2digits(val):
    if val < 10:
        return get1digit(val)
    if val < 20:
        return teens[val]
    else:
        return tens[val/10] + get1digit(val % 10)

def get3digits(val):
    n = val / 100
    mod = val % 100
    text =""
    if n > 0:
        text = get1digit(n) + "Hundred"
    if mod > 0:
        text += get2digits(mod)
    return text


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "parameters:"
        print "[1] - input file with numbers"
        quit()

    logging.basicConfig(stream=sys.stdout)
    logging.root.setLevel(logging.ERROR)
    fn = sys.argv[1]

    f = open(fn, "r")
    for line in f:
        if line == "":
            continue
        print number_to_text(int(line))
    f.close()

    exit(0)



