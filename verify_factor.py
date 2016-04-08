#!/usr/bin/env python3

import sys

if len(sys.argv) != 3:
    print("Run me as {} NUMBER FACTOR".format(sys.argv[0]), file=sys.stderr)
    sys.exit(2)

num, factor = int(sys.argv[1], base=0), int(sys.argv[2], base=0)
factor2 = num // factor
num2 = factor * factor2
if num == num2:
#    print("OK")
    sys.exit(0)
else:
#    print("NOT OK: {} {}".format(num, factor))
    sys.exit(1)
