#!/usr/bin/env python

import csv

# Convenience method to strip lines starting with a hashmark from a CSV file being read (assumed to be comments).
def decomment(csvfile):
    for row in csvfile:
        raw = row.split('#')[0].strip()
        if raw: yield raw

