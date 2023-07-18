#!/usr/bin/env python3

import sys, json, os

pipeline_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

config = json.load(open(os.path.join(pipeline_path, "config.json"), "r"))
obj = config
for arg in sys.argv[1:len(sys.argv)]:
	obj = obj[arg]
for key, val in obj.items():
	print(f"{key}=\"{val}\"")
