#!/usr/bin/env python

"""
Data structures for processing MC output.
"""

import collections

# Container for pseudoparticle point data
PPoint = collections.namedtuple("PPoint", ["r", "th", "phi", "e", "t"])

# Container for beginning and end points of trajectory (ie, a minimal trajectory)
MinTrajectory = collections.namedtuple("MinTrajectory", ["first", "last"])

