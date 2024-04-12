#!/usr/bin/env python3

import sys, re

t = sys.stdin.read().strip()
sci_pattern = '\d+\.\d+e-\d+'
print(re.sub(sci_pattern, lambda m: format(float(m.group()), '.10f'), t))
