#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 16:52:01 2022

@author: nico
"""
import re
fname = "220604_NR.tex"
with open(fname,"r") as f:
    text = f.read()

p = re.compile(r"\\nr{")

it = p.finditer(text)

comments = []
for match in it:
    cnt = 1
    for i, c in enumerate(text[match.start()+5:]):
        if c == "{":
            cnt += 1
        if c == "}":
            cnt -= 1
        if cnt == 0:
            break
    comments.append([match.start(), match.start() + 6 + i])